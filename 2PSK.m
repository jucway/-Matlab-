clear;
clc;
%输入信号
N = 1000; %发送序列长度Number of symbols
x = randi([0,1],1,N);
%生成[0,1]之间均匀分布的伪随机整数，作为1s时间内发送序列，发送速率为Rb=1000
code = 2*x - 1;   %产生双极性不归零码

fb = 1000; %发送端符号速率
fs = 16000; %滤波器采样频率
oversamp = fs / fb; %过采样率
delay = 5; %滤波器时延
alpha = 0.25; %滚降系数

%发送端脉冲成型

%平方根升余弦滤波器h_sqrt
h_sqrt = rcosine(1, oversamp, 'fir/sqrt', alpha, delay);   
sendsignal_oversample = kron(code, [1, zeros(1, oversamp-1)]); 
%首先对发送信号进行过采样
sendshaped = conv(sendsignal_oversample, h_sqrt);  %频域乘积，时域卷积
%sendshaped即为经过平方根升余弦脉冲成型后的时域波形
figure(1);
%画出平方根升余弦脉冲成型信号的时域波形
subplot(2,1,1),plot(sendshaped);axis([0 800 -1 1]);
title('发送端平方根升余弦脉冲成型后的信号时域波形');

f = (0:N-1)*fs/N;
x1 = fft(sendshaped, N); %信号的傅里叶变换
m1 = abs(x1); %取振幅
subplot(2,1,2),plot(f,m1);   %画出平方根升余弦脉冲成型信号的频域波形
title('发送端平方根升余弦脉冲成型后的频域波形');
xlabel('频率/Hz');
grid on;

%调制
fc = 4000; %载波频率
nt = 0:length(sendshaped) - 1;
carrier_wave1 = cos(2*pi*fc*nt/fs); %生成载波
modem_wave = sendshaped .* carrier_wave1; %调制（模拟相乘法）
figure(2);
subplot(2,1,1),plot(modem_wave);axis([0 800 -1 1]);    %画出已调信号的时域波形
title('发送端进行2ASK调制后的信号时域波形');

x2 = fft(modem_wave, N); %进行傅里叶变换
m2 = abs(x2); %取振幅
subplot(2,1,2),plot(f,m2); %画出已调信号的频域波形
xlabel('频率/Hz');
title('发送端进行2ASK调制后的信号频谱');
grid on;

%信道传输 加入高斯白噪声
snr = 10; %信噪比
TransSignal = awgn(modem_wave, snr,'measured','db'); %TransSignal为经过信道传输后得到的信号


%接收端相干解调
wave1 = TransSignal .* carrier_wave1; %乘以与发送端同频同相的载波 
figure(3);
subplot(2,1,1),plot(wave1); axis([0 800 -1 1]);  %画出接收端经过相乘器后的信号的时域波形
title('接收端进行解调后的信号时域波形');

x3 = fft(wave1, N); %信号的傅里叶变换
m3 = abs(x3); %取振幅
subplot(2,1,2),plot(f,m3); %画出接收端经过相乘器后的信号的频域波形
xlabel('频率/Hz');
title('接收端进行解调后的信号频谱');


%接收端低通滤波器
RecMatched1 = conv(wave1, h_sqrt);  %频域乘积，时域卷积
figure(4);
subplot(2,1,1),plot(RecMatched1);axis([0 800 -1 1]);%画出接收端接收后的信号的时域波形
title('接收端进行匹配滤波后的时域波形');

x4 = fft(RecMatched1, N); %信号的傅里叶变换
m4= abs(x4); %取振幅
subplot(2,1,2),plot(f,m4);   %画出接收端接收后的信号的频域波形
xlabel('频率/Hz');
title('接收端进行匹配滤波后的频谱');
grid on;

%接收端同步
SynPosi = delay * oversamp * 2 + 1;
SymPosi = SynPosi + (0:oversamp:(N-1) * oversamp);
RecSignal1 = RecMatched1(SymPosi);   %不发生倒π现象


%接收端采样 
RecBit1 = zeros(N); %初始化一个长度为 N 的全零数组，存放接收端抽样判决得到的单极性序列
for i = 1:N
    if RecSignal1(i)>0
        RecBit1(i)=1;
    elseif RecSignal1(i)<0
        RecBit1(i)=-1;
    end
end

figure(5);
subplot(2,1,1),stem(code,'.');axis([0 100 -1 1]);
title('发送端发送的序列');
subplot(2,1,2),stem(RecBit1,'.');axis([0 100 -1 1]);
title('接收端端抽样判决后得到的序列');
grid on;

%误码率
k = 0;
for i = 1:N
    if RecBit1(i) ~= code(i)
        k = k + 1;
    end
end
error1 = k / N;
sprintf('误码率：%2.2f%',error1);



%发生倒π情况

%接收端相干解调
%接收端相干载波是cos(wt+π)的情况
carrier_wave2 = cos(2*pi*fc*nt/fs + pi); %接收端在倒pi现象
wave2 = TransSignal .* carrier_wave2; %乘以与发送端同频同相的载波 
figure(6);
subplot(2,1,1),plot(wave2); axis([0 800 -1 1]);  %画出接收端经过相乘器后的信号的时域波形
title('接收端进行解调后倒π的信号时域波形');
 
x32 = fft(wave2, N); %信号的傅里叶变换
m32 = abs(x32); %取振幅
subplot(2,1,2),plot(f,m32);  %画出接收端经过相乘器后的信号的频域波形
xlabel('频率/Hz');
title('接收端进行解调后倒π的信号频谱');
grid on;

%接收端低通滤波器
RecMatched2 = conv(wave2, h_sqrt);  %频域乘积，时域卷积
figure(7);
subplot(2,1,1),plot(RecMatched2);axis([0 800 -1 1]);%画出接收端接收后的信号的时域波形
title('接收端进行匹配滤波后倒π的时域波形');
 
x42 = fft(RecMatched2, N); %信号的傅里叶变换
m42 = abs(x42); %取振幅
subplot(2,1,2),plot(f,m42);    %画出接收端接收后的信号的频域波形
xlabel('频率/Hz');
title('接收端进行匹配滤波后倒π的频谱');
grid on;

%接收端同步
RecSignal2 = RecMatched2(SymPosi);   %发生倒π现象

%接收端采样
RecBit2 = zeros(N); %初始化一个长度为 N 的全零数组，存放接收端抽样判决得到的单极性序列
for i = 1:N
    if RecSignal2(i)>0
        RecBit2(i)=1;
    elseif RecSignal2(i)<0
        RecBit2(i)=-1;
    end
end

figure(8);
subplot(2,1,1),stem(code,'.');axis([0 100 -1 1]);
title('发送端发送的序列');
subplot(2,1,2),stem(RecBit2,'.');axis([0 100 -1 1]);
title('接收端发送到π现象抽样判决后得到的序列');
grid on;

%误码率
j = 0;
for m = 1:N
    if RecBit2(m) ~= code(m)
        j = j + 1;
    end
end
error2 = j / N;
sprintf('误码率：%2.2f%',error2);
 



