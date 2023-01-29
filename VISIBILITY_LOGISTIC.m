%  LOAD THELFF DATA WITH ONE OR TWO FEEDBACKS
%
% CALCULATE THE PROBABILITY THAT GIVEN AN OP THE NEXT EVENT IS AN EXTREME
% EVENT.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
rmdir('output_data','s')
mkdir('output_data')
%%%%%%%%%
disp(' ')
disp('Andres Aragoneses and Nhat Nguyen code')
disp('----------------------')
disp('VISIBILITY LOGISTIC')
disp(date)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% GENERATE THE TIME SERIES %%%%
disp('----------------------')
disp('Generating the time series.');
%%%%%%%%%%%% 
L = 3000;    % Length of the time series
rmin = 3;
rmax = 4;
deltar = 0.0001; 
number = (rmax-rmin)/deltar+1;
skip = 500;   % Number of initial points removed from time series
rho = 0.23;
beta = 0.002;
K= 0.04;

for i=1:1:number;
    r(i) = rmin+(i-1)*deltar;

    %%%  Initial conditions
    ini = 0;
    fin = 1;
    x = (fin-ini).*rand(1) + ini;
%%%%%%%%%%%%%%%%%%%%%%           THE MAP        %%%%%%%%%%%%%%%%%%%%%
for j=1:1:L;
    x(j+1) = r(i)*x(j)*(1-x(j));         % LOGISTIC  % r from 3 to 4, but the window of periodicity starts from 3.7
%     x(j+1) = r(i)*min(x(j),1-x(j));       % TENT      % r from 1 to 2
%     x(j+1) = 1-r(i).*(sqrt(abs(x(j))));      % CUSP      % r from 0 to 1
%     x(j+1) = r(i)*sin(pi*x(j));        % SINE      % r from 0 to 1
%     x(j+1) = r(i)*x(j).*exp(-x(j));    % RICKERS   % r from 0 to 20
    % x(j+1) = rho + (K/(2*pi))*(sin(2*pi*x(j))+alpha*sin(4*pi*x(j)))+(beta*(1+1/K)*randn(1));
    % x(j+1) = rho + (K/(2*pi)) * ( sin(2*pi*x(j)) + alpha(i)*sin(4*pi*x(j))) + (beta*randn(1));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x(skip:length(x));

save('output_data/series.txt','x', '-ascii', '-append')

end

save('output_data/alpha.txt','r','-ascii','-append')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% CALCULATING THE WORDS %%%%

file = load('output_data/series.txt');
file2 = load('output_data/alpha.txt');

n = 1 ;    %% LAG TIME

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------------')
disp('Calculating the words.');

for index2 = 1:1:length(file2);
    parameter = file2(index2);
    ISI = file(index2,:);

clear a;                       % Clear any previous values
clear w;
clear A;
clear Prob;

%%%%%%%%%%%%%%%%%%%%%%%% .  CREATE THE WORDS OUT OF THE ISIs .   %%%%%%%%%%

w012=0;
w021=0;
w102=0;
w120=0;
w201=0;
w210=0;
regular = 0;
% T1 = 0;
% T2 = 0;

for i=1:1:length(ISI)-2;
    if ISI(i)<ISI(i+1) && ISI(i+1) < ISI(i+2);
        w012 = w012 +1;
        
    else
         if ISI(i)<ISI(i+2) && ISI(i+2) < ISI(i+1);
        w021 = w021 + 1;
         else
              if ISI(i+1)<ISI(i) && ISI(i) < ISI(i+2);
         w102 = w102 +1;
              else
                   if ISI(i+1)<ISI(i+2) && ISI(i+2) < ISI(i);
         w120 = w120 + 1;
                   else
                        if ISI(i+2)<ISI(i) && ISI(i) < ISI(i+1);
         w201 = w201 +1;
                        else
                             if ISI(i+2)<ISI(i+1) && ISI(i+1) < ISI(i);
         w210 = w210 + 1;
                             else
                                regular = regular + 1; % disp('equal IDIs');
                             end
                        end
                   end
              end
         end
    end
end

total = w012+w021+w102+w120+w201+w210;
one = w012/total;
two = w021/total;
three = w102/total;
four = w120/total;
five = w201/total;
six = w210/total;
% reg = regular;
seven = one + two + three;
eight = four + five + six;


words = [one two three four five six];
words1 = words;
words1(words1==0)=[];
entropy = -sum(words1.*log(words1)/log(6));
words2 = [two four one five three six];
Fo = 1;
FIM = 1;    % Fo*sum(diff(sqrt(words)).^2);
FIM_1 = 1;
FIM_2 = 1;
for i=1:1:6
    FIM_1 = Fo*sum(diff(sqrt(words2)).^2);
    if i >= 2
        FIM_2 = 1/2*sum(diff(sqrt(words2)).^2);        
    end
    FIM = FIM_1+ FIM_2;
end

P_R1 = (one + six)/2;          % 012, 210 
P_R2 = (two + three)/2;        % 021, 102 
P_R3 = (four + five)/2;        % 120, 201 
P_M1 = (one + six)/2;          % 012, 210
P_M2 = (four + two)/2;         % 120, 021
P_M3 = (five + three)/2;       % 201, 102
P_all = (one + two + three + four + five + six)/6; 

RV = ((one-P_R1).^2 + (six-P_R1).^2 + (five-P_R3).^2 + (four-P_R3).^2 + (two-P_R2).^2 + (three-P_R2).^2)/6;                        %Rotvar_psi1
RH = ((P_R1-P_all).^2 + (P_R2-P_all).^2 + (P_R3-P_all).^2)/2;                                                                      %Rothierarchy_psi2

MV = ((one-P_M1).^2 + (six-P_M1).^2 + (four-P_M2).^2)/6 + ((two-P_M2).^2 + (five-P_M3).^2 + (three-P_M3).^2)/6;                    %Mirrorvar_sigma1
MH = ((P_M1-P_all).^2 + (P_M2-P_all).^2 + (P_M3-P_all).^2)/3;                                                                      %Mirrorhierarchy_sigma2
  

P1 = [one-1/6 two+three-1/3];
PS1 = std(P1);

P2 = [six-1/6 four+five-1/3];
PS2 = std(P2);

ABS1 = abs((one-1/6))+abs((two-1/6))+abs((three-1/6));
ABS2 = abs((six-1/6))+abs((four-1/6))+abs((five-1/6));

ABS3 = abs((two-three));
ABS4 = abs((four-five));

%%% new defintions %%%
% V1 =1- abs((abs(one-1/6)-(abs(two-1/6+three-1/6)))/(abs(one-1/6)+abs(two-1/6+three-1/6)));%- abs((abs(six-1/6)-(abs(four-1/6)+abs(five-1/6)))/(abs(six-1/6)+abs(four-1/6)+abs(five-1/6)));      %V_alpha
% V2 =1- abs((abs(six-1/6)-(abs(four-1/6+five-1/6)))/(abs(six-1/6)+abs(four-1/6+five-1/6)));      %V_alpha

V1 = 1-(abs((abs(one-1/6)-(abs(two-1/6+three-1/6)))));%.*PS1;%- abs((abs(six-1/6)-(abs(four-1/6)+abs(five-1/6)))/(abs(six-1/6)+abs(four-1/6)+abs(five-1/6)));      %V_alpha
V2 = 1-(abs((abs(six-1/6)-(abs(four-1/6+five-1/6)))));%.*PS2;      %V_alpha

V23 = 1-abs((two-three))-abs((four-five));
V45 = 1-abs(abs(one-six))-abs(two-four)-abs(three-five);

% V23 = abs(two-three+four-five);%   abs(abs((two-three))+abs((four-five)));
% V45 = abs(abs(one-six)+abs(two-four)+abs(three-five));

%%% BANDT%%
Bbeta = abs(one-six);
Btau = abs(Bbeta-1/3);
Bgamma = abs(three+four- two-five);
Bdelta = abs(two+three-four-five);

 
%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE THE GRAY ZONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 1/6;
z = sqrt(m*(1-m)/total);
z1 = m+3*z;
z2 = m-3*z;

%%%%%%%%%%%%%%%%%%%%%

A = [parameter one two three four five six total entropy z1 z2 FIM regular];
B = [PS1 PS2 ABS1 ABS2 ABS3 ABS4 V23 V45 V1 V2 Bbeta Btau Bgamma Bdelta RV RH MV MH];
save('output_data/words.txt','A','-ascii','-append')
save('output_data/variance.txt','B','-ascii','-append')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------------')
disp('Plotting.');

file = load('output_data/words.txt');
file2 = load('output_data/variance.txt'); 
% current = 1:1:17;
% current = [10 20 30 40 50 60 70 80 90 100 110 120 130 140 150];
w1 = file(:,2);
w2 = file(:,3);
w3 = file(:,4);
w4 = file(:,5);
w5 = file(:,6);
w6 = file(:,7);
events = file(:,8);
PE = file(:,9);
error1 = file(:,10);
error2 = file(:,11);
fim = file(:,12);
reg = file(:,13);
current = file(:,1);
ps1 = file2(:,1);
ps2 = file2(:,2);
abs1 = file2(:,3);
abs2 = file2(:,4);
abs3 = file2(:,5);
abs4 = file2(:,6);
v23 = file2(:,7);
v45 = file2(:,8);
v1 = file2(:,9);
v2 = file2(:,10);
Bandt1 = file2(:,11);
Bandt2 = file2(:,12);
Bandt3 = file2(:,13);
Bandt4 = file2(:,14);
rv = file2(:,15);
rh = file2(:,16);
mv = file2(:,17);
mh = file2(:,18);

S1 = PE(1:length(PE)-1);
S2 = PE(2:length(PE));
dPE = S2-S1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% figures %%%

% figure(11)
% 
% subplot(221)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% l1 = plot(current,w1, 'b-','Markerfacecolor','b','linewidth',1);
% l2 = plot(current,w2, 'k-','Markerfacecolor','k','linewidth',1);
% l3 = plot(current,w3, 'g-','Markerfacecolor','g','linewidth',1);
% l4 = plot(current,w4, 'm-','Markerfacecolor','m','linewidth',1);
% l5 = plot(current,w5, 'c-','Markerfacecolor','c','linewidth',1);
% l6 = plot(current,w6, 'r-','Markerfacecolor','r','linewidth',1);
% e1 = plot(current,error1, 'k:','linewidth',1);
% e2 = plot(current,error2, 'k:','linewidth',1);
% legend('012','021','102','120','201','210')
% grid on
% ylabel('Words probabilities')
% xlabel('r')
% title('Logistics Map')
% xlim([rmin rmax])
% 
% subplot(222)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(current, v1,'Color','b','Markerfacecolor','b','linewidth',1);
% plot(current, v2,'Color','r','Markerfacecolor','r','linewidth',1);
% plot(current, v23, 'Color','m','Markerfacecolor','m','linewidth',1);
% plot(current, v45,'Color','k','Markerfacecolor','k','linewidth',1);
% title('Visibility')
% ylabel('Visibility')
% xlabel('r')
% legend('V_\alpha','V_\beta','V_{\delta}', 'V_{\rho}')
% grid on
% xlim([rmin rmax])
% 
% subplot(223)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(current, rv, 'Color','b','linewidth',1);
% plot(current, rh,'Color','r','linewidth',1);
% plot(current, mv,'Color','m','linewidth',1);
% plot(current, mh,'Color','c','linewidth',1);
% ylabel('Symmetry')
% xlabel('r')
% title('Rotational and Mirror Symmetry')
% legend('RV','RH','MV','MH')
% grid on
% xlim([rmin rmax])
% 
% subplot(224)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% % plot(current, v1,'-o','linewidth',1);
% % plot(current, v2,'-o','linewidth',1);
% % plot(current, v23,'-o','linewidth',1);
% % plot(current, v45,'-o','linewidth',1);
% plot(current, Bandt1,'linewidth',1);
% plot(current, Bandt2,'linewidth',1);
% plot(current, Bandt3,'linewidth',1);
% plot(current, Bandt4,'linewidth',1);
% title('Bandt')
% xlabel('r')
% legend('Interpreter','latex')
% legend('\beta', '\tau', '\gamma', '\delta')
% grid on
% xlim([rmin rmax])
% 
% figure(2)
% 
% subplot(2,1,1)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(PE, fim,'b.','markerfacecolor','b','linewidth',1);
% grid on
% ylabel('FIM')
% xlabel('PE')
% title('Logistics Map')
% 
% subplot(2,1,2)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(S1, dPE,'b.','markerfacecolor','b','linewidth',1);
% grid on
% ylabel('diff PE')
% xlabel('PE')
% title('Logistics Map')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%   LYAPUNOV EXPONENT FOR THE LOGISTIC MAP
% % 
% 
% lyap=zeros(1,1000);
% j=0;
% for(rr=rmin:deltar:rmax)
%     xn1=rand(1);
%     lyp=0;
%     j=j+1;
%     for(i=1:10000)
%         xn=xn1;
%         %logistic map
%         xn1=rr*xn*(1-xn);
%        %wait for transient
%        if(i>300)
%            % calculate teh sum of logaritm
%            lyp=lyp+log(abs(rr-2*rr*xn1));
%        end
%     end
%     %calculate lyapun
%     lyp=lyp/10000;
%     lyap(j)=lyp;
%     control(j)=rr;
% end
% 
% figure(3)
% 
% subplot(2,1,1)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(current,PE,'b-','markerfacecolor','b','linewidth',1);
% grid on
% ylabel('PE')
% xlabel('r')
% title('Logistics Map')
% xlim([rmin rmax])
% 
% subplot(2,1,2)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(control,lyap,'b-','markerfacecolor','b','linewidth',1);
% plot([rmin,rmax],[0,0], 'r-');
% grid on
% ylabel('\lambda')
% xlabel('r')
% title('Logistics Map')
% xlim([rmin rmax])


% figure(4)
% subplot(2,2,1)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(PE,w1, 'b.','Markerfacecolor','b','linewidth',1);
% plot(PE,w2, 'k.','Markerfacecolor','k','linewidth',1);
% plot(PE,w3, 'g.','Markerfacecolor','g','linewidth',1);
% plot(PE,w4, 'm.','Markerfacecolor','m','linewidth',1);
% plot(PE,w5, 'c.','Markerfacecolor','c','linewidth',1);
% plot(PE,w6, 'r.','Markerfacecolor','r','linewidth',1);
% ylabel('Probability')
% legend('w1','w2','w3','w4','w5','w6')
% 
% subplot(222)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(PE, v1,'b.','Markerfacecolor','b','linewidth',1);
% plot(PE, v2,'r.','Markerfacecolor','r','linewidth',1);
% plot(PE, v23,'m.','Markerfacecolor','m','linewidth',1);
% plot(PE, v45,'k.','Markerfacecolor','k','linewidth',1);
% legend('V_\alpha','V_\beta','V_{\delta}', 'V_{\rho}')
% grid on
% 
% subplot(2,2,3)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(PE, rv,'b.','linewidth',1);
% plot(PE, rh,'r.','linewidth',1);
% plot(PE, mv,'m.','linewidth',1);
% plot(PE, mh,'c.','linewidth',1);
% ylabel('Symmetry')
% xlabel('PE')
% legend('RV','RH','MV','MH')
% grid on
% 
% subplot(2,2,4)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(PE, Bandt1,'b.','linewidth',1);
% plot(PE, Bandt2,'r.','linewidth',1);
% plot(PE, Bandt3,'k.','linewidth',1);
% plot(PE, Bandt4,'m.','linewidth',1);
% xlabel('PE')
% legend('Interpreter','latex')
% legend('\beta', '\tau', '\gamma', '\delta')
% grid on

% figure(2)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(PE, v1,'b.','Markerfacecolor','b','linewidth',1);
% plot(PE, v2,'r.','Markerfacecolor','r','linewidth',1);
% plot(PE, v23,'m.','Markerfacecolor','m','linewidth',1);
% plot(PE, v45,'k.','Markerfacecolor','k','linewidth',1);
% legend('V_\alpha','V_\beta','V_{\delta}', 'V_{\rho}')
% grid on

% figure(5)
% 
% subplot(2,1,1)
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% l1 = plot(current,w1, 'b-o','Markerfacecolor','b','linewidth',1);
% l2 = plot(current,w2, 'k-o','Markerfacecolor','k','linewidth',1);
% l3 = plot(current,w3, 'g-o','Markerfacecolor','g','linewidth',1);
% l4 = plot(current,w4, 'm-o','Markerfacecolor','m','linewidth',1);
% l5 = plot(current,w5, 'c-o','Markerfacecolor','c','linewidth',1);
% l6 = plot(current,w6, 'r-o','Markerfacecolor','r','linewidth',1);
% e1 = plot(current,error1, 'k:','linewidth',1);
% e2 = plot(current,error2, 'k:','linewidth',1);
% legend('012','021','102','120','201','210')
% grid on
% ylabel('Probabilities')
% 
% subplot(2,1,2)
% hold all
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(current, v1,'-o','Color','b','Markerfacecolor','b','linewidth',1);
% plot(current, v2,'-o','Color','r','Markerfacecolor','r','linewidth',1);
% plot(current, v23,'-o', 'Color','c','Markerfacecolor','c','linewidth',1);
% plot(current, v45,'-o','Color','k','Markerfacecolor','k','linewidth',1);
% ylabel('Visibility')
% xlabel('Modulation amplitude')
% legend('V_\alpha','V_\beta','V_{\delta}', 'V_{\rho}')
% grid on

% subplot(3,1,3)
% hold all
% hold all
% set(gca,'Fontsize',24,'linewidth',1); box on
% plot(current, v1-v2,'-o','Color','b','Markerfacecolor','b','linewidth',0.5);
% % plot(current, v1+v2,'-o','Color','r','Markerfacecolor','r','linewidth',1);
% % plot(current, v23,'-o', 'Color','c','Markerfacecolor','c','linewidth',1);
% % plot(current, v45,'-o','Color','k','Markerfacecolor','k','linewidth',1);
% ylabel('Visibility')
% xlabel('Modulation amplitude')
% legend('V_\alpha - V_\beta')
% grid on


