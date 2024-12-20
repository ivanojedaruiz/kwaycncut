function FI=fft_filt_2(V,FB,sf);
% FI=fft_filt_2(V,FB,sf);
% fft-based filtering
% requires image to be called "V"
% and filters to be in FB
% sf is the subsampling factor =1
%
% FI is the result
% Jianbo Shi, 1997

[M1,M2,N3]=size(FB);
% prepare FFT of image for filtering
[N1,N2]=size(V);
I=zeros(size(V)+[M1-1 M2-1]);
I(1:N1,1:N2)=V;
N1s=length(1:sf:N1);
N2s=length(1:sf:N2);
IF=fft2(I);  %fast fourier on original image
FI=zeros(N1s,N2s,N3);

% apply filters
for n=1:N3
   %f=rot90(FB(:,:,n),2); %Rotate filter bank 180 degrees
   fF=fft2(FB(:,:,n),N1+M1-1,N2+M2-1);  %fast fourier on filter
   IfF=IF.*fF;  %multiply both fourier transforms
   If=real(ifft2(IfF));  %inverse fourier (real part?)
   If=If(ceil((M1+1)/2):ceil((M1+1)/2)+N1-1,ceil((M2+1)/2):ceil((M2+1)/2)+N2-1); 
   FI(:,:,n)=If(1:sf:N1,1:sf:N2);  %save in n dimensional matrix array
end

