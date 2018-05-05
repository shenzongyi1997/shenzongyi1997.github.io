%LS_ESPRIT ALOGRITHM
%DOA ESTIMATION BY LS_ESPRIT ALGORITHM
clear all;
close all;
%%�����еı��������
source_number=3;%��Դ��
sensor_number=8;%�����е���Ԫ��
m=7;%����Ԫ��
N_x=1024; %�źų���
snapshot_number=N_x;%������
w=[pi/15 pi/6 pi/3].';%�ź�Ƶ��
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2)+(2*pi*3e8)/w(3))/3;%�źŲ���
d=0.5*l;%��Ԫ���
snr=0;%�����/dB
source_doa=[59 45 5];%�źŵ�����Ƕ�
A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(1)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(2)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(3)*pi/180)/l)].';%��������
%A=(exp(-j*(0:sensor_number-1)*d*2*pi/l))'*sin(source_doa(1:source_number)*pi/180);
s=10.^(snr/20)*exp(j*w*[0:N_x-1]);%�����ź�
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%���˸�˹������������н����ź�





%%��ʱ��ʼESPRIT�ˣ�
x1=x(1:m,:);%����1���ܵ�����ʸ��
x2=x(2:(m+1),:);%����2���ܵ�����ʸ��
%�����������ģ�ͽ��кϲ�
X=[x1;x2];
%R=X*X'/snapshot_number;
R=X*X';
%��R��������ֵ�ֽ�
[U,S,V]=svd(R);
%R=R-S(2*m,2*m)*eye(2*m);%�������
%[U,S,V]=svd(R);
Us=U(:,1:source_number);%ȡ������Է�йص�������������
%disp('Us');
%disp(Us);
Us1=Us(1:m,:);
Us2=Us(m+1:2*m,:);
%���չ�ʽ�õ���ת�������M
psi=pinv(Us1)*pinv(Us1')*Us1'*Us2;%���ڷ���A�����з���B��ʹ�ã�A��B=B��A=I�����BΪA�������
%disp('M');
%disp(M);
%�Եõ�����ת���������������ֽ�
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

