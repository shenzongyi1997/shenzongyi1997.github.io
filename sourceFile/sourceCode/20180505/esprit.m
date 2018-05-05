%LS_ESPRIT ALOGRITHM
%DOA ESTIMATION BY LS_ESPRIT ALGORITHM
clear all;
close all;
%%把所有的变量都设好
source_number=3;%信源数
sensor_number=8;%总阵列的真元数
m=7;%子阵元数
N_x=1024; %信号长度
snapshot_number=N_x;%快拍数
w=[pi/15 pi/6 pi/3].';%信号频率
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2)+(2*pi*3e8)/w(3))/3;%信号波长
d=0.5*l;%阵元间距
snr=0;%信噪比/dB
source_doa=[59 45 5];%信号的入射角度
A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(1)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(2)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(3)*pi/180)/l)].';%阵列流型
%A=(exp(-j*(0:sensor_number-1)*d*2*pi/l))'*sin(source_doa(1:source_number)*pi/180);
s=10.^(snr/20)*exp(j*w*[0:N_x-1]);%仿真信号
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%加了高斯白噪声后的阵列接收信号





%%是时候开始ESPRIT了！
x1=x(1:m,:);%子阵1接受的数据矢量
x2=x(2:(m+1),:);%子阵2接受的数据矢量
%对两个子阵的模型进行合并
X=[x1;x2];
%R=X*X'/snapshot_number;
R=X*X';
%对R进行奇异值分解
[U,S,V]=svd(R);
%R=R-S(2*m,2*m)*eye(2*m);%降噪过程
%[U,S,V]=svd(R);
Us=U(:,1:source_number);%取出和鑫苑有关的两个特征向量
%disp('Us');
%disp(Us);
Us1=Us(1:m,:);
Us2=Us(m+1:2*m,:);
%按照公式得到旋转不变矩阵M
psi=pinv(Us1)*pinv(Us1')*Us1'*Us2;%对于方阵A，若有方阵B，使得：A·B=B·A=I，则称B为A的逆矩阵。
%disp('M');
%disp(M);
%对得到的旋转不变矩阵进行特征分解
[Vm,Dm]=eig(psi);
%disp('Dm');
%disp(Dm);
Dm=(diag(Dm)).';
%  theta = 0:1:180;
% Err = zeros(length(theta));
% idx1 = 0:m-1;
% idx2 = m:2*m-1:
% A = [idx1*d*exp(-2*pi*j/
Dmangle = -angle(Dm)/pi;
for theta1=1:90
    Err1(theta1) = (sin(theta1*pi/180)-Dmangle(1)).^(-2);
end
figure
 
    [doa,index] = max(Err1);
    Err1=Err1/doa;   
    plot(Err1);
     hold on
    plot([source_doa(1) source_doa(1)],[0 1]);
   for theta2=1:90
    Err2(theta2) = (sin(theta2*pi/180)-Dmangle(2)).^(-2);
   end
   [doa,index] = max(Err2);
   Err2 = Err2/doa;
   figure
    plot(Err2);
    hold on
    plot([source_doa(2) source_doa(2)],[0 1]);
for theta3=1:90
    Err3(theta3) = (sin(theta3*pi/180)-Dmangle(3)).^(-2);
end
figure

[doa,index] = max(Err3);
Err3=Err3/doa;
plot(Err3);
     hold on
    plot([source_doa(3) source_doa(3)],[0 1]);
doa=-asin(angle(Dm)/pi)*180/pi;
%disp('doa');
disp(doa);

