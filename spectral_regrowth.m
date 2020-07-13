function [x,y] = spectral_regrowth(M,b)
    phi1=pi/6;
    phi2=pi/4.5;
    h1=exp(-1j*pi*(0:M-1)*sin(phi1));
    h2=exp(-1j*pi*(0:M-1)*sin(phi2));
    H=[h1' h2']';
    P=H'*inv(H*H');
    m_order=4;
    b_in_s1=randi([0 1], log2(m_order)*1000,1);
    b_in_s2=randi([0 1], log2(m_order)*1000,1);
    s1 = qammod(b_in_s1, m_order, 'InputType', 'bit', 'UnitAveragePower', true);
    s2 = zeros(1000,1);
    s1_upsampled=upsample(s1,4);
    filtr=rcosdesign(0.5,8,4);
    s1_tilde=filter(filtr,1,s1_upsampled);
    s2_upsampled=upsample(s2,4);
    s2_tilde=filter(filtr,1,s2_upsampled);
    s_tilde=[s1_tilde s2_tilde]';
    x=P*s_tilde;
    range=linspace(-1/M,1/M,2^b);
    y_real=[];
    y_imag=[];
    y=zeros(M,4000);
    for n=1:size(x,2)
        y_real=[];
        y_imag=[];
        for m=1:size(x,1)
            realpart=real(x(m,n));
            imagpart=imag(x(m,n));
            [valreal,realidx]=min(abs(range-realpart));
            [valimag,imagidx]=min(abs(range-imagpart));
            y_real=[y_real range(realidx)];
            y_imag=[y_imag range(imagidx)];
        end
        y(:,n)=complex(y_real,y_imag);
    end
end