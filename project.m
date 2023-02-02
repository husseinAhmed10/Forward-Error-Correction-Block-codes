clc
clear
close all
x=[0 1 0 1 1 0 1 0 1 0 1 0];
n=31;
SNR=3;

% %without FEC
% mod= BPSK_Modulation(x);
% y_awgn = awgn(mod,10);
% z=Receiver_without_FEC(y_awgn)
% demod=BPSK_Demodulation(z);

% %without FEC ->> PAM
%mod= fn_PAM_Modulation(x);
%y_awgn = awgn(mod,5);
%demod=fn_Pam_demodulation(y_awgn);

% %with FEC %Rept ->> BPSK
% y_rep= Repetition_Encoder(n,x);
% mod = BPSK_Modulation(y_rep);
% y_awgn =  awgn(mod,1);
% demod = BPSK_Demodulation(y_awgn);
% z=Repetition_Decoder(n,demod);

% %with FEC %Rept ->>PAM
% y_rep= Repetition_Encoder(n,x);
% mod = fn_PAM_Modulation(y_rep);
% y_awgn =  awgn(mod,1);
% demod = fn_Pam_demodulation(y_awgn);
% z=Repetition_Decoder(n,demod);

% %with FEC %Hamming ->> BPSK
% [y_rep,padding]= Hamming_Encoder(n,x);
% mod = BPSK_Modulation(y_rep);
% y_awgn =  awgn(mod,5);
% demod = BPSK_Demodulation(y_awgn);
% z=Hamming_Decoder(n,demod,padding);

% %with FEC %Hamming ->>PAM
% [y_rep,padding]= Hamming_Encoder(n,x);
% mod = fn_PAM_Modulation(y_rep);
% y_awgn =  awgn(mod,10);
% demod = fn_Pam_demodulation(y_awgn);
% z=Hamming_Decoder(n,demod,padding);

%%BER 
%x=BER_without_FEC ("BPSK", [0 1 0 1 1 1 0 0], 10, 0);
%y=BER_with_FEC ("Repition code","BPSK", [0 1 0 1 1 1 0 0], 10, 0, 3);
%z=BER_with_FEC ("Hamming","BPSK", [0 1 0 1 1 1 0 0], 10, 0, 7);



%AWGN (Didn't use and used awgn) (i.e. both random generator)
        function [Signal_Dis]=AWGN_BPSK_OR_PAM(SNR_DB,EncodedSignal)
        % Signal_Dis= awgn(EncodedSignal,SNR_DB) ; % or can use awgn
        noise=(randn(1,length(EncodedSignal))); %% generating a random Bits (this is the Noise)
        N0=1/(10^(SNR_DB/10)); %%this is the N0 
        sigma=sqrt(N0/2); %%this is the sigma of the Noise 
        noise=noise*sigma; %%Then muliply the sigma by the Noise I have
        Signal_Dis=(EncodedSignal+noise); %%Then addedd the Noise to the Encoded Bits
        
        end

%%Reciever of BPSK without AWGN
        function[DecodedSignal]= Receiver_without_FEC(Bits)
        
            for i=1:length(Bits)
                if(Bits(i)<0)
                    DecodedSignal(i)=-1;
                elseif(Bits(i)>=0)
                    DecodedSignal(i)=1;
                end
            end
            
            %Demodulation BPSK
%             for i=1:length(DecodedSignal)
%                 if(DecodedSignal(i)==-1)
%                     DecodedSignal(i)=0;
%                 end  
%             end
        end
        
        %k=4;
        %Type = 1; %(1-> rep, 2-> ham)
%%functions used in PAM
function [PAM_Decoder1]=fn_decimal_binary(PAM_Decoder)
%Pam from decimal to binary%%%%%%%%%%%%%%%%%%%
Length_Pam_decoder=length(PAM_Decoder);
PAM_Decoder1=[];   
   p=1;
    for j=1:1:Length_Pam_decoder
        if(PAM_Decoder(j)==0)
            PAM_Decoder1(p)=0;
            PAM_Decoder1(p+1)=0;
            PAM_Decoder1(p+2)=0;
       elseif(PAM_Decoder(j)==1)
           PAM_Decoder1(p)=0;
           
            PAM_Decoder1(p+1)=0;
           
            PAM_Decoder1(p+2)=1;
            
       elseif(PAM_Decoder(j)==2)
             PAM_Decoder1(p)=0;
            PAM_Decoder1(p+1)=1;
            PAM_Decoder1(p+2)=0;
            
       elseif(PAM_Decoder(j)==3)
             PAM_Decoder1(p)=0;
            PAM_Decoder1(p+1)=1;
            PAM_Decoder1(p+2)=1;
       elseif(PAM_Decoder(j)==4)
             PAM_Decoder1(p)=1;
            PAM_Decoder1(p+1)=0;
            PAM_Decoder1(p+2)=0;
       elseif(PAM_Decoder(j)==5)
             PAM_Decoder1(p)=1;
            PAM_Decoder1(p+1)=0;
            PAM_Decoder1(p+2)=1;
       elseif(PAM_Decoder(j)==6)
             PAM_Decoder1(p)=1;
            PAM_Decoder1(p+1)=1;
            PAM_Decoder1(p+2)=0;
       elseif(PAM_Decoder(j)==7)
          PAM_Decoder1(p)=1;
            PAM_Decoder1(p+1)=1;
            PAM_Decoder1(p+2)=1;
       end
        p=p+3;
    end
end

function [message]=fn_conversion_Bin_decimal(Message)
m=length(Message);
message=[];
a=mod(m,3);
 if(a==0)
    p=1;
    for j=1:1:m/3
        message(j)=(1*Message(p+2))+(2*Message(p+1))+(4*Message(p));
        p=p+3;
    end
   else
    disp("There are single or double bit in array are remaining");   
end
end
    

   %%Fair 
 function [snr_fair] = Fair(Type,Snr,k,n)

        snr_fair = snr;

            if (Type==1) %Repition
                %EncodedSignal_out=EncodedSignal./(sqrt(n));
                snr_fair = snr_fair -10*log(n);
            elseif (Type==2) % hamming
              %EncodedSignal_out=EncodedSignal ./(sqrt(n/k));
              snr_fair = snr_fair -10*log(n/k) ;
            end
            
        end
     
        
   %%BPSK Modulation and Demodulation
        
        function [Modulated]=BPSK_Modulation(stream)
        
            Modulated = stream;
            for i=1:length(stream)
                if(stream(i)==0)
                    Modulated(i)=-1;
                end
            end
        
        end
      
        function [Demodulated]=BPSK_Demodulation(Bits)
        
        Demodulated = Bits;
        for i=1:length(Bits)
    
            if(Bits(i)<0)
        
                Demodulated(i)=0;
    
            elseif(Bits(i)>=0)
        
                Demodulated(i)=1;
   
            end

        end
        end
        
   %%PAM Modulation and Demodulation
   function [mess3]=fn_Pam_demodulation(Noisy)
Length_of_received_signal=length(Noisy);
mess3=[];
   for j=1:1:Length_of_received_signal
       if(Noisy(j)==-7 ||Noisy(j)<-7 || (Noisy(j) > -7 && Noisy(j) < -6) )
            mess3(j)=0;
       elseif(Noisy(j)==-5 || (Noisy(j) > -6 && Noisy(j) < -4))
            mess3(j)=1;
       elseif(Noisy(j)==-1 || (Noisy(j) > -2 && Noisy(j) < 0))
             mess3(j)=2;
       elseif(Noisy(j)==-3 || (Noisy(j) > -4 && Noisy(j) < -2))
             mess3(j)=3;
       elseif(Noisy(j)==7 ||Noisy(j)>6)
             mess3(j)=4;
       elseif(Noisy(j)==5 || (Noisy(j) > 4 && Noisy(j) < 6))
             mess3(j)=5;
       elseif(Noisy(j)==1 || (Noisy(j) > 0 && Noisy(j) < 2))
             mess3(j)=6;
       elseif(Noisy(j)==3 || (Noisy(j) > 2 && Noisy(j) < 4))
           mess3(j)=7;
       end
       
   end
end

function [mess1]=fn_PAM_Modulation(message)
p=length(message);
mess1=[];
   for j=1:1:p
       if(message(j)==0)
            mess1(j)=-7;
       elseif(message(j)==1)
            mess1(j)=-5;
       elseif(message(j)==2)
             mess1(j)=-1;
       elseif(message(j)==3)
             mess1(j)=-3;
       elseif(message(j)==4)
             mess1(j)=7;
       elseif(message(j)==5)
             mess1(j)=5;
       elseif(message(j)==6)
             mess1(j)=1;
       elseif(message(j)==7)
           mess1(j)=3;
       end
       
   end
end
        
        
 %%Repition
        
function [Encoded_Repetition] = Repetition_Encoder (n,message)

Encoded_Repetition = zeros(1,n*length(message)); % empty array

k=1;

for i=1:n: (length(Encoded_Repetition)-1)
    
    for j=1:1:n
        Encoded_Repetition(1,i+j-1) = message(k);
    end
    k = k+1 ;
    
end

%Encoded_BPSK_Repetition = real( pskmod(Encoded_Repetition,2,pi) );

end

function [Decoded_BPSK] = Repetition_Decoder (n,Recieved)


Decoded_BPSK = zeros(1,length(Recieved)/n);


for i=1:(length(Recieved)/n)
    
   num_ones = 0;
   num_zeros = 0;
    for j=1:1:n 
        
        if Recieved(n*(i-1)+j)== 1
            num_ones = num_ones +1;
        else
            num_zeros = num_zeros +1;
        end
    end     
        if(num_ones>num_zeros)
            Decoded_BPSK(1,i)=1;
        elseif (num_ones<num_zeros)
            Decoded_BPSK(1,i)=0;
        else
             Decoded_BPSK(1,i)=1;   %n is even and equal number of errors for both which is n/2 so majority rule can't be achieved
        end         
    
end

end

%%Hamming

function [EncodedMessage,PaddingZero] = Hamming_Encoder (n,message)

 %Encoded_BPSK_Hamming 1*n
            %message 1*k
            %n = {7, 15, 31}
            %k = {4, 11, 26}
            %m = {3, 4 , 5 }
            %Generating matrix is k*n
            if n==7
                m=3;
                GenMat= gen2par(hammgen(m));
                k=n-m;
            elseif n==15
                m=4;
                GenMat= gen2par(hammgen(m));
                k=n-m;
            elseif n==31
                m=5;
                GenMat= gen2par(hammgen(m));
                k=n-m;
            end
            
            lengthOfMessage=length(message);
            %%no of padding zeros
            PaddingZero = mod((k-mod(lengthOfMessage,k)),k);
            padded_array = [zeros(1,PaddingZero),message];
            length_new_message = PaddingZero+lengthOfMessage;   %length of new encoded message without hamming
            EncodedMessage=zeros([1 round(length_new_message*(n/k))]); %number of i/ps included
            
            encoded_diff_size= rem(transpose(fliplr(reshape(padded_array,k,[])))*GenMat,2);
            EncodedMessage=reshape(encoded_diff_size.',1,[]);
            
            % %Another method
            % for i = 1:(length_new_message/k)
            %   EncodedMessage(n*(i-1)+1:n*(i-1)+n)=rem(padded_array(k*(i-1)+1:k*(i-1)+k)*GenMat,2);
            % end


end

function [DecodedMessage] = Hamming_Decoder (n,Recieved,padding)

%Encoded_BPSK_Hamming 1*n
            %message 1*k
            %n = {7, 15, 31}
            %k = {4, 11, 26}
            %m = {3, 4 , 5 }
            %Generating matrix is k*n
            
            % dmin= 3 -> correct 1 bit
            % all syndrome can correct 1 bit and detect 2 bit, but cannot correct them
            %so more than 1 bit error will be a repeated syndrome so that we take the
            %the less number of errors which is one bit
            %number of states of syndrome equals 2^n-k = 2^m
            
            if n == 7
                GenMat= gen2par(hammgen(3)); %%%%%%%hammgen(m)
                %H_trans = [eye(3);GenMat(:,1:3)];
                H_trans = transpose(hammgen(3));
                %Error=[zeros(1,7);eye(7)];
                %syndrome_table = Error*H_trans;
                syndrome_table=syndtable(transpose(H_trans));
                k=4;
                m=3;
            elseif n == 15
                GenMat = gen2par(hammgen(4));%%%%%%%hammgen(m)
                %H_trans = [eye(3);GenMat(:,1:4)];
                H_trans = transpose(hammgen(4));
                %Error=[zeros(1,15);eye(15)];
                %syndrome_table = Error*H_trans;
                syndrome_table=syndtable(transpose(H_trans));
                k=11;
                m=4;
            elseif n == 31
                GenMat = gen2par(hammgen(5));%%%%%%%hammgen(m)
                %H_trans = [eye(3);GenMat(:,1:5)];
                H_trans = transpose(hammgen(5));
                %Error=[zeros(1,31);eye(31)];
                %syndrome_table = Error*H_trans;
                syndrome_table=syndtable(transpose(H_trans));
                k=26;
                m=5;
            end
            
            
            % syndrome = mod(Recieved*H_trans,2);
            % for i=1:1: (n+1)
            % if syndrome_table(i,(1:m))== syndrome
            %     %disp("test1")
            %     error_found =  Error(i,(1:n));
            % end
            % end
            % corrected = mod(Recieved + error_found, 2);%+ mod(syndrome*((H_trans)^(-1)),2) ;
            % Decoded_BPSK_Hamming = corrected(k:n) ;
            
            %The above code for correcting at the same time if you don't want to
            %calculate BER
            
            
        %   DecodedMessage=zeros([1,length(Recieved)*(k/n)]);

            if n ==15
   
                k1 = 5;

            elseif n==7
    
                k1=4;
            elseif n==31
   
                k1=6;
            end
            for i= 1:(length(Recieved)/n)
                  VHT = bi2de(rem(flip(transpose(fliplr(reshape(Recieved,n,[]))))*H_trans,2),'left-msb'); % Convert to decimal.
                  e = syndrome_table(1+VHT,:); 
                  CorrectedMessage=rem(e+flip(transpose(fliplr(reshape(Recieved,n,[])))),2);    
   
            end
            CorrectedMessage=CorrectedMessage(:,[m+1:end]);
            CorrectedMessage=flip(CorrectedMessage);
            DecodedMessage=reshape(CorrectedMessage.',1,[]);
            DecodedMessage = DecodedMessage(1,padding+1:end);
            
end

%%BER

function [BER]= BER_without_FEC (BPSK_PAM, bit_stream, SNR_max, SNR_min)

  
if BPSK_PAM == "BPSK"
                    
                   
    Stream_modulated = BPSK_Modulation(bit_stream);    
                    
                
else
    
    Stream_modulated = fn_PAM_Modulation(bit_stream);
                    
               
end
                counter=1;
                Error_Bit=0;
                BER=[];
                for j=SNR_min:((SNR_max-SNR_min)/10):SNR_max
                    
                    BER(counter)=0;
                    Signal_Dis = awgn(Stream_modulated,j);
                    %counter = counter + 1;
                    % test=Receiver_without_FEC(app,Signal_Dis);
                    if BPSK_PAM == "BPSK"
                        
                        Received_Signal = BPSK_Demodulation(Signal_Dis);
                        
                    else
                        Received_Signal=fn_Pam_demodulation(Signal_Dis);
                    end
                    
                    
                    
                    for i=1:1:length(Received_Signal)
                        
                        if Received_Signal(i)~=bit_stream(i)
                            Error_Bit=Error_Bit+1;
                            
                        end
                        
                    end
                    Error_Bit =Error_Bit/length(Received_Signal);
                    BER(counter) = Error_Bit;
                    counter=counter+1;
                    
                end
                
                
                P=0.5*erfc(sqrt(10.^((SNR_min:1:SNR_max)/10)));
                
                semilogy((SNR_min:((SNR_max-SNR_min)/10):SNR_max),BER,'LineWidth',2);
              %  legend(app.UIAxes,{'without FEC'},'Location','southoutside');
           
end

function [BER] = BER_with_FEC (type,BPSK_PAM, bit_stream, SNR_max, SNR_min, n)

if type == "Repition code"  %repition
                    
                    Encoded_bits = Repetition_Encoder(n,bit_stream);
                    
                    if BPSK_PAM == "BPSK"
                        
                        Stream_modulated = BPSK_Modulation(Encoded_bits);
                        
                    else
                        Stream_modulated = fn_PAM_Modulation(Encoded_bits);
                        % test_hw(app);
                        %    uialert(uifigure,'File not found','Invalid File');
                    end
                    
                    counter=1;
                    Error_Bit=0;
                    BER=[];
                    for j=SNR_min:((SNR_max-SNR_min)/10):SNR_max
                        
                        BER(counter)=0;
                        Signal_Dis = awgn(Stream_modulated,j);
                        %counter = counter + 1;
                        % test=Receiver_without_FEC(app,Signal_Dis);
                        if BPSK_PAM == "BPSK"
                            
                            Received_Signal = BPSK_Demodulation(Signal_Dis);
                            
                            
                        else
                            
                            Received_Signal = fn_Pam_demodulation(Signal_Dis);
                            
                            
                        end
                        
                        Received_Signal=Repetition_Decoder(n,Received_Signal);
                        
                        
                        
                        for i=1:1:length(Received_Signal)
                            
                            if Received_Signal(i)~=bit_stream(i)
                                Error_Bit=Error_Bit+1;
                                
                            end
                            
                        end
                        Error_Bit =Error_Bit/length(bit_stream);
                        BER(counter) = Error_Bit;
                        counter=counter+1;
                        
                    end
                    
                    
                    P=0.5*erfc(sqrt(10.^(SNR_min:((SNR_max-SNR_min)/10):SNR_max)));
                    
                    semilogy((SNR_min:((SNR_max-SNR_min)/10):SNR_max),BER,'LineWidth',2);       %Plot the line code
                 %   legend(app.UIAxes,{'Repition'},'Location','southoutside');
                
                else    %hamming
                    
                    if n ~= 7 && n ~= 15 && n ~= 31
                        disp("re enter n of hamming");
                        return;
                    end
                    
                   
                    [Encoded_bits,padding] = Hamming_Encoder(n,bit_stream);
                    
                    if BPSK_PAM == "BPSK"
                        
                        Stream_modulated = BPSK_Modulation(Encoded_bits);
                        
                    else
                        Stream_modulated = fn_PAM_Modulation(Encoded_bits);
                      
                    end
                    
                  counter=1;
                  Error_Bit=0;
                  BER=[];
                    
                    
                    for j=SNR_min:((SNR_max-SNR_min)/10):SNR_max
                        
                        BER(counter)=0;
                        Signal_Dis = awgn(Stream_modulated,j);
                        %counter = counter + 1;
                        % test=Receiver_without_FEC(app,Signal_Dis);
                        if BPSK_PAM == "BPSK"
                            
                            Received_Signal = BPSK_Demodulation(Signal_Dis);
                            
                            
                        else
                            
                            Received_Signal = fn_Pam_demodulation(Signal_Dis);
                            
                            
                        end
                        
                        Received_Signal=Hamming_Decoder(n,Received_Signal,padding);
                        
                        
                        
                        for i=1:1:length(Received_Signal)
                            
                            if Received_Signal(i)~=bit_stream(i)
                                Error_Bit=Error_Bit+1;
                                
                            end
                            
                        end
                        Error_Bit =Error_Bit/length(Received_Signal);
                        BER(counter) = Error_Bit;
                        counter=counter+1;
                        
                    end
                    
                    
                    semilogy((SNR_min:((SNR_max-SNR_min)/10):SNR_max),BER,'LineWidth',2);       %Plot the line code
                 %   legend(app.UIAxes,{'Hamming'},'Location','southoutside');
                    
                end

end



















