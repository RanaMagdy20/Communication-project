clear;
%%%%%%part1
M =  input('Enter input matrix(each column represents a signal): ');
[n m] = size(M);
T = input ('Enter the max time duration: ')
Tb=T/n;
figure
  for i=(1:m)
  t_M= transpose(M(:,i));
   t=0:Tb:T;
  amp= [t_M 0];
  subplot(m,1,i)
  stairs(t,amp,'b','LineWidth',2)
  title_text1 = sprintf('Input signal s_%d(t)' , i);
  title_text2 = sprintf('s_%d(t)' , i);
  title(title_text1);
  xlabel('time(sec)');
  ylabel(title_text2);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%energy of each symbols using the constellation diagram.%%%%%%%%%%%%%%%%%%
E=zeros(1,m);
E(1)=coef(1,1)^2;
Energy=sprintf("Energy of S1=%d",E(1))
 for i=2:m %where i is the signal number
for j=1:b %where j is the basis number
   E(i)=E(i)+coef(i,j)^2;
end
 energy=sprintf("Energy of S%d=%d",i,E(i))
end
%%%%%%
%part2
PNRZamp = [1 -1; 1 -1];
Tb = 1;
[phi,coef] =Gram_schmidt(PNRZamp,Tb);
constellation_dig(1,2,coef);

%%%%Generate stream of random bits
random_bits= randi([0 1],1,100000);
bitrate=1;
n = 6; %how much sample point will be appeared for each bit.
T=length(random_bits)/bitrate; %  total time needed for the bit stream.
N= n*length(random_bits);  %total sample point needed for whole bit stream
ts=T/N;       %time for each sample point
t=0:ts:T-ts;
Tb=1;    % Bit rate


%%%% For Polar non return to zero

polar_non_return_to_zero = zeros(1,length(t));
for i=1:length(random_bits)
  if random_bits(i)==1
   polar_non_return_to_zero((i-1)*n+1:i*n) = 1;
  else  polar_non_return_to_zero((i-1)*n+1:i*n) = -1;
  end
end

Tx=polar_non_return_to_zero ;

%%%%making wgn

rxdec=zeros(1,length(random_bits));
 BER = zeros(1,9); %Size allocating the BER matrix
for k= 1:9
s=0;
Tb=1;
SNR= [-10, -8, -6, -4 ,-2, 0, 2, 4, 6]; %in dB
Eb= 10*log(Tb); %Eb in dB (Eb=0 db)
Nodb= Eb - SNR(k);     %energy of noise in dB
No = 10^ (Nodb/10);    %power of noise in watt

noise_real= wgn(1, N , No/2 , 'linear');      %"linear" to state that N0/2 is in watts not db
noise_img= 1j*wgn(1, N , No/2 , 'linear');
noise= noise_real+noise_img;
Rx = Tx + noise; %adding the wgn to the transmitted signal

%%%%%drawing constellation diagram with noise

PNRZamp = [1 -1; 1 -1];
[Phi,coef]= Gram_schmidt (PNRZamp,1);
 figure 
constellation_dig(1,2,coef);
hold on
  scatter(real(Rx), imag(Rx), 'b.');
hold off
  title(['Constellation Digram at Eb/No = ' num2str(SNR(k)) 'dB']);

axis ([-10 10 -10 10]);
%%%%%decoding recived data, each bit is made from 6 points

%%desicion maker
c=1;
for u=3:6:length(Rx)

    if(Rx(u)>0)
            rxdec(c)=1;
    else
            rxdec(c)=0;
    end
   c=c+1;
end

% BER by comparing
for L=1:length(random_bits)
    if(rxdec(L)~=random_bits(L))      % compare the recieved decoded data to the original bits transmitted
    s=s+1;
    end
end
BER(k)=s/length(random_bits);   %Probability of error of the decoded signal by comparing


end

figure()

 semilogy(SNR,BER, 'color', 'b','linewidth',3);
title('The Probability of Error Due to Noise');
ylabel('BER');
xlabel('Eb/No (dB)');
hold on;

%%Theoretically BER
x = sqrt(10.^(SNR/10));
BERT = 0.5* erfc(x);
semilogy(SNR,BERT, 'r-*' ,'linewidth',2);
legend('BER by Comparing','Theortical BER');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%constellation diagram%%%%%%%%%%%%%%%%%%%
function constellation_dig(b,m,coef) % b=no. of basis function & m=no. of signals
S=zeros(m,b);
if (b == 1) %in case of 1 basis function
    for i=1:m
        S(i,1)=coef(i,1);
       % plot(S(i,1),0,"o") %ploting 1D graph for phi1 on x axis
        hold on
    end
         title("Constellaton Digram"); %adjusting title
         xlabel("\phi_1"); %adjusting x label
         ylabel("\phi_2"); %adjusting y label
         xlim([min(min(coef))-1 max(max(coef))+1]); %adjusting X-axis
         ylim([min(min(coef))-1 max(max(coef))+1]); %adjusting Y-axis
         grid on;
     elseif(b==2)
        for i=1:m
          for j=1:b
          S(i,j)=coef(i,j);
          end
           % plot(S(i,1),S(i,2),"o")
            hold on
        end
         title("Constellaton Digram"); %adjusting title
         xlabel("\phi_"+int2str(1)); %adjusting x label
         ylabel("\phi_"+int2str(2)); %adjusting y label
         xlim([min(min(coef))-1 max(max(coef))+1]); %adjusting X-axis
         ylim([min(min(coef))-1 max(max(coef))+1]); %adjusting Y-axis
         grid on;
        figure(4)
             elseif (b==3) %in case of 3 basis function
             for i=1:m
                for j=1:b
                  S(i,j)=coef(i,j);
                end
               % plot3(S(i,1),S(i,2),S(i,3),"x");
                hold on;
                end
                title("Constellaton Digram"); %adjusting title for the graph
                xlabel("\phi_"+int2str(1));
                ylabel("\phi_"+int2str(2));
                zlabel("\phi_"+int2str(3));
                xlim([min(min(coef))-1 max(max(coef))+1]); %adjusting X-axis
                ylim([min(min(coef))-1 max(max(coef))+1]); %adjusting Y-axis
                grid on;
                else
                  warning("can't plot more than 3 basis fanctions in constellation digram");

end

end
function [phi,coef] = Gram_schmidt(M,T)
[n m] = size(M);
Tb=T/n;
m2=M.^2;
m2_sum=sum(m2,1);
Energy = m2_sum.*Tb;
Energy_s=sqrt(Energy);
phi=zeros(n,m); %intiallizing phi matrix to store results
coef= zeros(m,m);
phi(:,1)=M(:,1)/Energy_s(1); %phi=signal/root of energy
coef(1,1)=Energy_s(1);

    for j= 2 : m
        g = M(:,j);

          for k = 1:j-1 % filling coef above s(j,j)
          coef(k,j) = (M(:,j)' * phi(:,k))*Tb;
          g = g - coef(k,j)*phi(:,k); %g=g-coef*phi
          end
            energy_g = sum(g.^2)*Tb;
        if (energy_g == 0 || energy_g < 1e-5)
            continue;
        end
            phi(:,j) = g/sqrt(energy_g);
    end
          for z = 2 : m %calculate the diagonal coefficients
            coef(z,z)= (M(:,z)' * phi(:,z))*Tb;
          end
phi=phi(:,any(phi)); %remove zeros in phi
coef2=coef';
coef=coef2(:,any(coef2)); %remove zeros in coef
[a b]=size(phi);
figure(2);
for c=(1:b)
  t_phi= transpose(phi(:,c));
  t=0:Tb:T;
  amp_P= [t_phi 0];
  %subplot(b,1,c)
  stairs(t,amp_P,'b','LineWidth',2)
  title_text1 = sprintf('basis signal s_%d(t)' , c);
  title_text2 = sprintf('phi_%d(t)' , c);
  title(title_text1);
  xlabel('time(sec)');
  ylabel(title_text2);
end
end