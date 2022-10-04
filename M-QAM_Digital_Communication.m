close all;
clear all;

blockLength = 1000;
nBlocks=10000;
M = 16;
sqM = sqrt(M);
EbdB = 1.0:1.0:12.0;
Eb = 10.^(EbdB/10);
n = log2(M);
Es = n*Eb;
No = 1;
SNR = Es/No;
SNRdB = 10*log10(SNR);
IdxI = randi([0,sqM-1],1,blockLength);
IdxQ = randi([0,sqM-1],1,blockLength);
Sym = (2*IdxI-(sqM-1))+1j*(2*IdxQ-(sqM-1));
SER = zeros(1,length(EbdB));

for blk = 1:nBlocks
    noise = sqrt(No/2)*(randn(1,blockLength)+1j*randn(1,blockLength));
    for K = 1:length(EbdB)
        A = sqrt(Es(K)*3/2/(M-1));
        TxSym = A*Sym;
        RxSym = TxSym + noise;
        EqSym = RxSym/A;
        DecIdxI = MQAM_DECODER(real(EqSym),M);
        DecIdxQ = MQAM_DECODER(imag(EqSym),M);
        SER(K) = SER(K)+sum(or(DecIdxI~=IdxI,DecIdxQ~=IdxQ));
    end
end


SER = SER/blockLength/nBlocks;

semilogy(SNRdB,SER,'g  s','linewidth',3.0,'MarkerFaceColor','g','MarkerSize',9.0);
hold on;
semilogy(SNRdB,4*(1-1/sqrt(M))*qfunc(sqrt(3*SNR/(M-1))),'b -  ','linewidth',3.0);
axis tight;
grid on;
legend('SER','SER Theory')
xlabel('SNR (dB)');
ylabel('SER');
title('SER vs SNR(dB) for Digital Comm with M-QAM');


