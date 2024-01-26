# 稳态大作业实现思路

## 潮流计算

### 牛拉法

要点是数据类型的转换和失配量的生成

```MATLAB
clc
clear;

mpc=case39;
node=length(mpc.bus(:,1));% number of nodes
nbranch=length(mpc.branch(:,1));% number of branch
tic
%% Generate Y matrix
yN=makeYbus(mpc);
G=real(yN);
B=imag(yN);
%% Process PQ, PV nodes
Vdelta=find(mpc.bus(:,2)==3);% Vdelta
PQ=find(mpc.bus(:,2)==1);   % PQ bus
PV=find(mpc.bus(:,2)==2);   % PV bus
n0=find(mpc.bus(:,2)<3);    

%% Get information from mpc.
nVdelta=length(Vdelta);      % number of vdelta node
nPQ=length(PQ);              % number of PQ nodes
nPV=length(PV);                 % number of PV nodes
V=mpc.bus(:,8);
delta=deg2rad(mpc.bus(:,9));  % angle(bus给的是角度，需要换算成弧度)
Pgen=zeros(node,1);
Qgen=zeros(node,1);
Pgen(mpc.gen(:,1))=mpc.gen(:,2);
Qgen(mpc.gen(:,1))=mpc.gen(:,3);
P_net=(Pgen-mpc.bus(:,3))./mpc.baseMVA; % 节点功率净输入量
Q_net=(Qgen-mpc.bus(:,4))./mpc.baseMVA;
tol_max=1e-8;    % tolerance of the solution

%% Iterative method Newton Raphson
tolerance=1; 
count=0;
start=datetime("now");
while (tolerance>tol_max)
    [P_cal,Q_cal]=cal_power(G,B,node,delta,V); % this function calculates the power
    [mismatch_PQ]=mismatch_Power(P_net,Q_net,P_cal,Q_cal,PQ,n0); % this function gets the mismatch power
    % Jacobian Matrix
    [J]=gen_Jacobi(nPQ,nPV,G,B,V,delta,PQ,n0,Vdelta); 
    mismatch_X=J\mismatch_PQ;
    ddelta=mismatch_X(1:node-1);
    dV=mismatch_X(node:end);
    %Update V and delta
    delta(n0)=delta(n0)+ddelta;
    V(PQ)=V(PQ)+dV;
    tolerance=max(abs(mismatch_PQ));
    count=count+1;
end
endl=datetime("now");
time=endl-start;
%% 计算线路功率和节点功率

% sdelta=-(delta*180/pi+3);
% delta=sdelta*pi/180;
sdelta=delta.*180/pi;
Vm = V.*exp(1j*delta);
I=yN*Vm;
S=Vm.*conj(I);
P_node=zeros(1,nPV+1);
Q_node=zeros(1,nPV+1);

% PV和Vdelta节点功率
for k=1:length(mpc.gen(:,1))
    P_node(k)=real(S(k))+mpc.gen(k,2);
    Q_node(k)=imag(S(k))+mpc.gen(k,3);
end
PQ_res=P_node+Q_node.*1j;
toc
disp("=========节点功率========")
disp(PQ_res.*mpc.baseMVA);
% 支路功率
for k=1:node
    for m=1:node
        Ss(k,m)=Vm(k)*conj((Vm(k)-Vm(m))*yN(k,m));
   end
end
[l1,u1]=lu(Ss);
disp("==========支路功率===========")
% Define active and reactive power flows
Pij=real(u1);
Qij=imag(u1);
Pji=real(l1);
Qji=imag(l1);
disp("============frombus to tobus===============")
disp(sparse(Pij+1j.*Qij).*mpc.baseMVA);
disp("============tobus to frombus===============")
disp(sparse(Pji+1j.*Qji).*mpc.baseMVA);


disp("============节点电压=============")
disp(sparse(V))
disp("============电压相角=============")
disp(sparse(sdelta))

disp("===========总损耗===========")
disp(-sum(Ss(:)).*mpc.baseMVA);
disp("============================")
disp(['运行时间: ',num2str(toc)]);
% 
% tic
% result1 = runpf(mpc,mpoption('pf.alg','NR','pf.tol',1e-8));
% toc
```
