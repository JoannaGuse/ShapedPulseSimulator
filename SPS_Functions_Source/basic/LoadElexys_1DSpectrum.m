function [spec]= LoadElexys_1DSpectrum(pathname, file)
%loads 1D spectrum data from .DTA and .DCS file.
%Syntax:
% [spec]= LoadElexys_1DSpectrum(pathname, file)
%
%_________outputs_________
% spec.fuw=microwave frequence
% spec.B0=center fiels
% spec.x= x axis-B field (Gauss)
% spec.y=spectrum

M=1;

file=strrep(file,'.DTA','');
file=strrep(file,'.dta','');
fid = fopen( [ pathname file '.DTA'], 'r','ieee-be.l64') ;% spc
%check existance of file
if(fid == -1)
    return;
end
%load DTA file and read into DSC file.
[Z,N]= fread(fid,inf,'float64');
fclose(fid);
fid = fopen([ pathname file '.DSC'], 'r');% par
[S,~]= fscanf(fid,'%c');
fclose(fid);

% %Parse string S, the delimiter is ascii 10, this is the token delimiter
Delim = setstr(10);
VB=[];
while(length(S) > 0 )
    [token,S] = strtok(S,Delim);
    VB = str2mat(VB,token);
end

GST = get_vb(VB,'XMIN');
GSI = get_vb(VB,'XWID');
N = get_vb(VB,'XPTS');
dx = GSI/(N-1);
x=GST + dx*(0:(N-1));
dim2=get_vb(VB,'YTYP');
if ~strncmp(dim2,'NODATA',3),
    ymin = get_vb(VB,'YMIN'); ywid = get_vb(VB,'YWID'); M = get_vb(VB,'YPTS');
    dy = ywid/(M-1);
    y=ymin + dy*[0:(M-1)];
end;

if length(Z)==2*N*M,
    z0=Z;
    Z=zeros(N,1);
    for k=1:N*M,
        Z(k)=z0(2*k-1)+1i*z0(2*k);
    end;
end;
Z=reshape(Z,N,M);
Z=permute(Z,[2,1]);
%% change x axis to MHz
B=(x*1e-4);
B0 =double(get_vb(VB,'A1CT'))*1e4;
fuw = double(get_vb(VB,'MWFQ'));

spec.fuw=fuw;
spec.B0=B0;
spec.x=B;
spec.y=Z;
end
%_______________________________________________
