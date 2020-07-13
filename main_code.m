%% Problem 1 Question 1
N=16;
phi=-pi/6;
lambda=1; 
d_set=[lambda/2, 2*lambda, lambda/8];
theta=-pi/2:0.01:pi/2;
gain = zeros(3,length(theta)); 
for d=1:length(d_set)
    gain = zeros(1,length(theta)); 
    w_BF = exp(-1j*2*pi*(0:N-1)'*d_set(d)*sin(phi)/lambda)/sqrt(N);
    for t = 1:length(theta)
        a_theta = exp(-1j*2*pi*(0:N-1)'*d_set(d)*sin(theta(t))/lambda)/sqrt(N);
        gain(t) = abs(w_BF'*a_theta)^2;
    end
    gainlog=10*log10(gain);
    polarplot(theta, gainlog, 'LineWidth', 1)
    hold on
end
legend('d=\lambda/2', 'd=2\lambda', 'd=\lambda/8', 'FontSize', 15)
rlim([-10 0])
thetalim([-90 90])
title('Normalized beam pattern', 'FontSize', 15)
%% Problem 2 Question 1
Nt=32;
dt=1/2;
dr=1/2;
L=4;
Pt=0.001;
B=200*10^3;
No=-174;
No=10^((No-30)/10);
% No=-204;
Nr=4:2:32;
c_r=zeros(length(Nr),1000);
c_s=zeros(length(Nr),1000);
c_r_fin=zeros(1,length(Nr));
c_s_fin=zeros(1,length(Nr));
for r=1:length(Nr)
    for m=1:1000
        Hs=zeros(Nr(r),Nt);
        Hr=randn(Nr(r),Nt);
        for l=1:L
            thetai=unifrnd(-pi/2,pi/2);
            phii=unifrnd(-pi/2,pi/2);
            alphai=normrnd(0,sqrt(max(Nt,Nr(r))));
            arx_thetai = exp(-1j*2*pi*(0:Nr(r)-1)'*dr*sin(alphai))/sqrt(Nr(r));
            atx_phii = exp(-1j*2*pi*(0:Nt-1)'*dt*sin(phii))/sqrt(Nt);
            Hs=Hs+alphai*arx_thetai*atx_phii';
        end
        kr=rank(Hr);
        ks=rank(Hs);
        [Us, Ss, Vs]=svd(Hs);
        [Ur, Sr, Vr]=svd(Hr);
        for k=1:kr
            c_r(r,m)=c_r(r,m)+log2(1+(Nr(r)*Nt*Pt*abs(Sr(k,k))^2/(kr*B*No)));
        end
        c_r(r,m)=B*c_r(r,m);
        for k=1:ks
            c_s(r,m)=c_s(r,m)+log2(1+(Nr(r)*Nt*Pt*abs(Ss(k,k))^2/(ks*B*No)));
        end
        c_s(r,m)=B*c_s(r,m);
    end
    c_r_fin(r)=mean(c_r(r,:));
    c_s_fin(r)=mean(c_s(r,:));
end
%% 
figure
semilogy(Nr,c_r_fin,Nr,c_s_fin, 'linewidth',2)
legend('Rich Scattering Capacity','Sparse Scattering Capacity'), title('\fontsize{14}Channel Capacity'), xlabel('\fontsize{12}Nr'), ylabel('\fontsize{12}Rate in bps')
%% Problem 2 Question 2
Nt=32;
dt=1/2;
dr=1/2;
L=4;
Pt=0.001;
B=200*10^3;
No=-174;
No=10^((No-30)/10);
% No=-204;
Nr=4:2:32;
achieveable_r=zeros(length(Nr),1000);
achieveable_s=zeros(length(Nr),1000);
achieveable_r_fin=zeros(1,length(Nr));
achieveable_s_fin=zeros(1,length(Nr));
achieveable_r_unk=zeros(length(Nr),1000);
achieveable_s_unk=zeros(length(Nr),1000);
achieveable_r_unk_fin=zeros(1,length(Nr));
achieveable_s_unk_fin=zeros(1,length(Nr));
ftr_angle=-pi/2:pi/Nt:-(-pi/2+pi/Nt);
ftr=exp(-1j*2*pi*(0:Nt-1)'*dt*sin(ftr_angle))/sqrt(Nt);
for r=1:length(Nr) 
    wtr_angle=-pi/2:pi/Nr(r):-(-pi/2+pi/Nr(r));
    wtr=exp(-1j*2*pi*(0:Nr(r)-1)'*dt*sin(wtr_angle))/sqrt(Nr(r));
    for m=1:1000
        Hs=zeros(Nr(r),Nt);
        Hr=randn(Nr(r),Nt);
        for l=1:L
            thetai=unifrnd(-pi/2,pi/2);
            phii=unifrnd(-pi/2,pi/2);
            alphai=normrnd(0,sqrt(max(Nt,Nr(r))));
            arx_thetai = exp(-1j*2*pi*(0:Nr(r)-1)'*dr*sin(alphai))/sqrt(Nr(r));
            atx_phii = exp(-1j*2*pi*(0:Nt-1)'*dt*sin(phii))/sqrt(Nt);
            Hs=Hs+alphai*arx_thetai*atx_phii';
        end
        kr=rank(Hr);
        ks=rank(Hs);
        [Us, Ss, Vs]=svd(Hs);
        [Ur, Sr, Vr]=svd(Hr);
        for k=1:2
            achieveable_r(r,m)=achieveable_r(r,m)+log2(1+(Nr(r)*Nt*Pt*abs(Sr(k,k))^2/(2*B*No)));
        end
        achieveable_r(r,m)=B*achieveable_r(r,m);
        for k=1:2
            achieveable_s(r,m)=achieveable_s(r,m)+log2(1+(Nr(r)*Nt*Pt*abs(Ss(k,k))^2/(2*B*No)));
        end
        achieveable_s(r,m)=B*achieveable_s(r,m);
        y_abs_r=[];
        y_abs_s=[];
        for i=1:size(ftr,2)
            for j=1:size(wtr,2)
                y_r=wtr(:,j)'*Hr*ftr(:,i);
                y_abs_r_temp=abs(y_r)^2;
                y_s=wtr(:,j)'*Hs*ftr(:,i);
                y_abs_s_temp=abs(y_s)^2;
                y_abs_r=[y_abs_r y_abs_r_temp];
                y_abs_s=[y_abs_s y_abs_s_temp];
            end
        end
        yrarray=sort(y_abs_r,'descend');
        ysarray=sort(y_abs_s,'descend');
        for k=1:2
            achieveable_r_unk(r,m)=achieveable_r_unk(r,m)+log2(1+(Nr(r)*Nt*Pt*yrarray(k)/(2*B*No)));
        end
        achieveable_r_unk(r,m)=B*achieveable_r_unk(r,m);
        for k=1:2
            achieveable_s_unk(r,m)=achieveable_s_unk(r,m)+log2(1+(Nr(r)*Nt*Pt*ysarray(k)/(2*B*No)));
        end
        achieveable_s_unk(r,m)=B*achieveable_s_unk(r,m);
    end
    achieveable_r_fin(r)=mean(achieveable_r(r,:));
    achieveable_s_fin(r)=mean(achieveable_s(r,:));
    achieveable_r_unk_fin(r)=mean(achieveable_r_unk(r,:));
    achieveable_s_unk_fin(r)=mean(achieveable_s_unk(r,:));
end
%% 
ytick1=10^7:5*10^7:5*10^8;
yticklab=({'0.6x10^8','1.1x10^8','1.6x10^8','2.1x10^8','2.6x10^8','3.1x10^8','3.6x10^8'});
figure
semilogy(Nr,c_r_fin, 'linewidth',2)
hold on
semilogy(Nr,c_s_fin,'linewidth',2)
hold on
semilogy(Nr,achieveable_r_fin,'linewidth',2)
hold on
semilogy(Nr,achieveable_s_fin,'linewidth',2)
hold on
semilogy(Nr,achieveable_r_unk_fin,'linewidth',2)
hold on
semilogy(Nr,achieveable_s_unk_fin,'linewidth',2)
% ylim('auto')
yticks(ytick1)
yticklabels(yticklab)
legend('Rich Scattering Capacity','Sparse Scattering Capacity','Rich Scattering Achievable (Ch known)','Sparse Scattering Achievable (Ch known)','Rich Scattering Achievable (Ch unknown)','Sparse Scattering Achievable (Ch unknown)', 'FontSize',12)
title('Capacity, Channel known/unknown achievable rate', 'FontSize', 12)
xlabel('\fontsize{12}Nr'), ylabel('\fontsize{12}Rate in bps')
%% Problem 3 question 1,2
M_values=[8,32];
b_values=[4,12];
colors=['b','r';'k','m'];
phi_range=-pi/2:0.01:pi/2;
for m=1:2
    for b=1:2
        [x,y]=spectral_regrowth(M_values(m),b_values(b));
        figure(1)
        [pxx,ww]=periodogram(y(1,:));
        plot(ww/pi,10*log10(pxx),colors(m,b), 'linewidth',1.5)
        hold on
        legend('M=8 b=4','M=8 b=12','M=32 b=4','M=32 b=12', 'FontSize',15)
        title('Impact of DAC quantization in the frequency domain','FontSize',15)
        xlabel('\fontsize{13}Normalized Frequency'), ylabel('\fontsize{13}Power Spectral Density')
        phi_range=-pi/2:0.01:pi/2;
        a=zeros(M_values(m),length(phi_range));
        for p=1:length(phi_range)
            a(:,p)=exp(-1j*pi*(0:M_values(m)-1)'*sin(phi_range(p)));
        end
        z=y;
%         clear g
        for p=1:length(phi_range)
            G=[];
            for n=1:4000
                G=[G abs(a(:,p)'*z(:,n))^2];
            end
            g(p)=rms(G);
        end
        figure(2)
        plot(phi_range*180/pi,10*log10(g), 'linewidth',2)
        hold on
        legend('M=8 b=4','M=8 b=12','M=32 b=4','M=32 b=12', 'FontSize',15)
        title('Impact of DAC quantization in the angular domain','FontSize',15)
        xlabel('\fontsize{13}Angle in degrees'), ylabel('\fontsize{13}Angular response')
    end
end
%% Problem 3 question 3,4
beta=[0,-133];
for m=1:2
    for be=1:2
        [x,y]=spectral_regrowth(M_values(m),b_values(b));
        z=zeros(M_values(m),4000);
        for n=1:4000
            z(:,n)=x(:,n)+beta(be)*x(:,n).*abs(x(:,n)).^2;
        end
        figure(3)
        [pxx,ww]=periodogram(z(1,:));
        plot(ww/pi,10*log10(pxx),colors(m,be))
        hold on
        legend('M=8 \beta_3=0','M=8 \beta_3=-133','M=32 \beta_3=0','M=32 \beta_3=-133', 'FontSize',15)
        xlabel('\fontsize{13}Normalized Frequency'), ylabel('\fontsize{13}Power Spectral Density')
        title('Impact of PA non-linearity in the frequency domain','FontSize',15)
        a=zeros(M_values(m),length(phi_range));
        for p=1:length(phi_range)
            a(:,p)=exp(-1j*pi*(0:M_values(m)-1)'*sin(phi_range(p)));
        end 
        for p=1:length(phi_range)
            G=[];
            for n=1:4000
                G=[G abs((a(:,p)'*z(:,n)))^2];
            end
            g(p)=rms(G);
        end
        figure(4)
        plot(phi_range*(180/pi),10*log10(g), 'linewidth',2)
        hold on
        legend('M=8 \beta_3=0','M=8 \beta_3=-133','M=32 \beta_3=0','M=32 \beta_3=-133', 'FontSize',15)
        title('Impact of PA non-linearity in the angular domain','FontSize',15)
        xlabel('\fontsize{13}Angle in degrees'), ylabel('\fontsize{13}Angular response')
    end
end