%% MRI DTMREIT analysis for the Pinapple
%% Ainslie Johnstone

cd ~/'OneDrive - University College London'/1_DTMREIT/
TE=[ 5, 7.32, 9.64, 11.92, 14.82, 16.6, 18.92, 21.24, 23.56, 25.88, 28.2, 30.52, 32.84, 35.16, 37.48, 39.8, 42.12, 44.44] ;
addpath(genpath('~/matlab/spm12'))
cd ~/'OneDrive - University College London'/1_DTMREIT/Pineapple/AP_stim/Mag
  AllMag=niftiread('AP_stim_4Dmag.nii.gz');
  info=niftiinfo('AP_stim_4Dmag.nii.gz');

  PosMag=AllMag(:,:,:,1:2:35);
  NegMag=AllMag(:,:,:,2:2:36);
  


%% The T2 star loop 
warning('off','all')
[sx,sy,sz,e]=size(PosMag);
for x=1:sx
    for y=1:sy
        for z=1:sz
            
                   thisvox=squeeze(PosMag(x,y,z,:));
   
            T2star= curvefit(TE, thisvox);
            T2sMap(x,y,z)=T2star;
        end
        disp(strcat( num2str(((y/sy*100)/sx)+(((x-1)/sx)*100)), '% finished'))
    end
end

T2sMap=single(T2sMap); 


  info.Filename=('/Users/ainsliej/OneDrive - University College London/1_DTMREIT/Pineapple/AP_stim/Mag/T2s_Pos.nii')
  info.Datatype='single';
  info.ImageSize=[84,96,52];
  niftiwrite(T2sMap, 'T2s_Pos', info);
    
%% Make a Bz image at each echo- LR version

% First make the complex image for the positive and negative current
% directions
% complex=mag*exp(1i*phase) %perhaps its got to be .multiplication 
% complex_dif=complex_pos.*conj(complex_neg)
% then for the RHS of equation angle(complex_diff) is equal to the arg part

% make a R2* image using an MPM analysis sequence, and add in the R2* map
% into the analysis for calculation of the weighting. 
% Perhaps spatially smooth the Bz image if there are any extreme values 

%clear all
cd ~/'OneDrive - University College London'/1_DTMREIT/Pineapple/LR_stim
addpath(genpath('~/matlab/spm12'))
TE=[5, 7.32, 9.64, 11.92, 14.82, 16.6, 18.92, 21.24, 23.56, 25.88, 28.2, 30.52, 32.84, 35.16, 37.48, 39.8, 42.12, 44.44] ;
Ctimes= (TE-3.5)/1000; % I am not sure on this, need to check when the current reached peak, but cant on Mac
 info=niftiinfo('../AP_stim/Fieldmaps/Bz_e1.nii');


for echo = 2:18
    if echo==1
         accrued_phase=[];
    end 
  phase=niftiread(strcat('Phase/LR_Pinapple_Phase_e',num2str(echo),'.nii.gz'));
  mag=niftiread(strcat('Mag/LR_Pinapple_Mag_e',num2str(echo),'.nii.gz'));
  phase = rescale(phase,-pi, pi,'InputMin', 0,'InputMax', 4096); % This is the correct range for Siemens phase images.
  phase_pos=squeeze(phase(:,:,:,1));
  phase_neg=squeeze(phase(:,:,:,2));
  mag_pos=double(squeeze(mag(:,:,:,1)));
  mag_neg=double(squeeze(mag(:,:,:,2)));
  
  complex_pos = mag_pos .* exp(1i * phase_pos);
  complex_neg = mag_neg .* exp(1i * phase_neg);
  complex_dif = complex_pos .* conj(complex_neg);
  
  accrued_phase = unwrap(cat(4, accrued_phase, complex_dif), [], 4);
    accrued_phase = accrued_phase(:,:,:,end);
    
    % Bz = [cycles during current] / [time current applied] / 2 (i.e. average) / Gamma
    Bz = accrued_phase / (2 * pi)  /      Ctimes(echo)      / 2                / 42.58e6; % Gamma = 42.58 MHz / T
    %Bz = 1/ (2 * 42.58e6 * Ctimes(echo)) * angle(complex_dif) / (2 * pi); % OLD VERSION
    %Gamma = 42.58 MHz / T
    % Bz = [cycles during current] / [time current applied] / 2 (i.e. average) / Gamma

    info.Filename=strcat('/Users/ainsliej/OneDrive - University College London/1_DTMREIT/Pineapple/LR_stim/Fieldmaps/Bz_e',num2str(echo),'.nii.gz');
    info.ImageSize=[84,96,52];
    info.PixelDimensions=[2.0833,2.0833,3.55];
    info.Datatype='double';
    name=strcat('Bz_e',num2str(echo));
    cd Fieldmaps
    niftiwrite(Bz,name ,info);
    clear Bz complex_dif complex_pos complex_neg phase mag phase_pos phase_neg mag_pos mag_neg
   
end 

%   Bz(1).dat=single(nansum(AllBz, 4)); 
%   Bz(1).mat=Bz.mat;
%   Bz(1).dim=Bz.dim;
%   info.Filename=('/Users/ainsliej/OneDrive - University College London/1_DTMREIT/Pineapple/AP_stim/AP_meanBz.nii.gz');
%   niftiwrite(meanBz(1).dat, 'AP_meanBz', info);   
%     
%% Make a Bz image at each echo- AP version

clear all
cd ~/'OneDrive - University College London'/1_DTMREIT/Pineapple/AP_stim
addpath(genpath('~/matlab/spm12'))
TE=[5, 7.32, 9.64, 11.92, 14.82, 16.6, 18.92, 21.24, 23.56, 25.88, 28.2, 30.52, 32.84, 35.16, 37.48, 39.8, 42.12, 44.44] ;
Ctimes= (TE-3.5)/1000; % I am not sure on this, need to check when the current reached peak, but cant on Mac
 info=niftiinfo('../AP_stim/Fieldmaps/Bz_e1.nii');
 accrued_phase=[];

for echo = 1:18
  phase=niftiread(strcat('Phase/AP_Pinapple_Phase_e',num2str(echo),'.nii.gz'));
  mag=niftiread(strcat('Mag/AP_Pinapple_Mag_e',num2str(echo),'.nii.gz'));
  phase = rescale(phase,-pi, pi,'InputMin', 0,'InputMax', 4096); % This is the correct range for Siemens phase images.
  phase_pos=squeeze(phase(:,:,:,1));
  phase_neg=squeeze(phase(:,:,:,2));
  mag_pos=double(squeeze(mag(:,:,:,1)));
  mag_neg=double(squeeze(mag(:,:,:,2)));
  
  complex_pos = mag_pos .* exp(1i * phase_pos);
  complex_neg = mag_neg .* exp(1i * phase_neg);
  complex_dif = complex_pos .* conj(complex_neg);
  
 % accrued_phase = unwrap(cat(4, accrued_phase, complex_dif), [], 4);
  %  accrued_phase = accrued_phase(:,:,:,end);
    
    % Bz = [cycles during current] / [time current applied] / 2 (i.e. average) / Gamma
    %Bz = accrued_phase / (2 * pi)  /      Ctimes(echo)      / 2                / 42.58e6; % Gamma = 42.58 MHz / T
    Bz = 1/ (2 * 42.58e6 * Ctimes(echo)) * angle(complex_dif) / (2 * pi); % OLD VERSION
    %Gamma = 42.58 MHz / T
    % Bz = [cycles during current] / [time current applied] / 2 (i.e. average) / Gamma

    info.Filename=strcat('/Users/ainsliej/OneDrive - University College London/1_DTMREIT/Pineapple/AP_stim/Fieldmaps/Bz_e',num2str(echo),'.nii.gz');
    info.ImageSize=[84,96,52];
    info.PixelDimensions=[2.0833,2.0833,3.55];
    info.Datatype='double';
    name=strcat('Bz_e',num2str(echo));
    cd Fieldmaps
    niftiwrite(Bz,name ,info);
    clear Bz complex_dif complex_pos complex_neg phase mag phase_pos phase_neg mag_pos mag_neg
   
end 


%% Combine echos with optimal weightings
cd ~/'OneDrive - University College London'/1_DTMREIT/Pineapple/LR_stim
addpath(genpath('~/matlab/spm12'))
T2s=niftiread('Mag/T2s_Pos.nii');
Mask=niftiread('Pinapple_mask_MREIT.nii.gz');
T2s=T2s.*Mask;
sz=size(T2s);
phi=zeros(sz(1), sz(2),sz(3),18);
sumphi=zeros(sz(1), sz(2),sz(3));
TE=[5, 7.32, 9.64, 11.92, 14.82, 16.6, 18.92, 21.24, 23.56, 25.88, 28.2, 30.52, 32.84, 35.16, 37.48, 39.8, 42.12, 44.44] ;
Ctimes= (TE-3.5); % I am not sure on this, need to check when the current reached peak, but cant on Mac
for e = 1:18
phi(:,:,:,e)=(Ctimes(e).*Ctimes(e)).*exp(-(2.*Ctimes(e))./(T2s));
sumphi=sum(phi, 4);
end


for echo = 1:18
    
    wf=phi(:,:,:,echo)./sumphi;
  Bz=spm_vol(strcat('Fieldmaps/Bz_e',num2str(echo),'_mas.nii.gz'));
  info=niftiinfo(char(strcat('Fieldmaps/Bz_e',num2str(echo),'_mas.nii.gz')));
  
  if echo==1 
      thisBz=Bz.dat;
      wBz=thisBz.*wf;
      AllBz=wBz;
  else 
      thisBz=Bz.dat;
      wBz=thisBz.*wf;
      AllBz=cat(4,AllBz,wBz);
  end
  
 
end

% AllBz=AllBz(:,:,:,[1:14,16:18]); %This is for the AP condition. where
% echo 15 is really weird!
  meanBz(1).dat=double(nansum(AllBz, 4)); 
  meanBz(1).mat=Bz.mat;
  meanBz(1).dim=Bz.dim;
  info.Filename=('/Users/ainsliej/OneDrive - University College London/1_DTMREIT/Pineapple/LR_stim/LR_meanBz.nii.gz');
  niftiwrite(meanBz(1).dat, 'LR_meanBz', info);  

  
  %% ROAST recipe 
  % I have been maually changing the eTMS file in the roast folder 
  % Also I have changed a bunch of things in electrode placement
  % Also final visualisation is broken, but it has run correctly 
  
  
 % rmpath(genpath('~/matlab/spm12'))
 %   cd ~/'OneDrive - University College London'/1_DTMREIT/Pineapple/AP_stim/roast
  addpath(genpath('/Users/ainsliej/matlab/roastV2.7.1'))
  cd /Users/ainsliej/matlab/roastV2.7.1
roast_PA('AProast/Pinapple_AP.nii', {'custom1', -1, 'custom2', 1}, ...
    'conductivities',struct('white',[1], 'gray',[1],'csf',[1], 'bone',[1], 'skin',[1]),...
    'electype',{'disc','disc'}, 'elecsize',{[10,2],[10,2]}, 'meshoptions',struct('radbound',4,'maxvol',8)) 
  %'zeropadding', 10,
  
    addpath(genpath('/Users/ainsliej/matlab/roastV2.7.1'))
  cd /Users/ainsliej/matlab/roastV2.7.1/LRroast
roast_PA('Pinapple_LR.nii', {'custom1', -1, 'custom2', 1}, ...
    'conductivities',struct('white',[1], 'gray',[1],'csf',[1], 'bone',[1], 'skin',[1]),...
    'electype',{'disc','disc'}, 'elecsize',{[10,2],[10,2]}, 'meshoptions',struct('radbound',4,'maxvol',8)) 



 %% Toolbox AP Jz
 
clear all 
cd /Users/ainsliej/'OneDrive - University College London'/1_DTMREIT/Pineapple/AP_stim/
addpath('~/OneDrive - University College London/1_DTMREIT/MRCI_toolbox')
VoxelSize = [3.55, 2.0833, 2.0833];
CurrAmpl = 0.002;       % injected current amplitude
Mask=int32(niftiread('Pinapple_mask_MREIT_RAS.nii'));
Bz=niftiread('AP_meanBz_RAS.nii.gz');
J0=niftiread('AP_J0image.nii');
% J0(:,:,:,1)=J0(:,:,:,1)/0.00355;
% J0(:,:,:,2)=J0(:,:,:,2)/0.0020833;
% J0(:,:,:,3)=J0(:,:,:,3)/0.0020833;
Bz(isnan(Bz))=0;
J0(isnan(J0))=0;
Mask(isnan(Mask))=0;
Mask=double(Mask);
Bz=double(Bz);
szJ0=size(J0);
szBz=size(Bz);
add=szBz(2)-szBz(1);
Bz=[Bz; zeros(add,szBz(2),szBz(3))];
J0=[J0; zeros(add,szJ0(2),szJ0(3), szJ0(4))];
Mask=[Mask; zeros(add,szBz(2),szBz(3))];


for thisSlice=1:szBz(3)
recon_parameters = reconstruction_parameters('projected current density type-1',...
    'VoxelSize',VoxelSize,'Mask',Mask(:,:,thisSlice));
thisJp=mrci_projected_current_density1(squeeze(J0(:,:,thisSlice,:)),squeeze(Bz(:,:,thisSlice)),recon_parameters) ;  
Jp(:,:,thisSlice,:)=thisJp;
disp(strcat('Just finished slice... ',num2str(thisSlice)))
end

info=niftiinfo('AP_J0image.nii');
info.Filename=('/Users/ainsliej/OneDrive - University College London/1_DTMREIT/Pineapple/AP_stim/AP_Jp.nii.gz');
info.ImageSize=[84,84,96,3];
niftiwrite(Jp, 'AP_Jp', info);


 
 %% Toolbox LR Jz
 
clear all 
cd /Users/ainsliej/'OneDrive - University College London'/1_DTMREIT/Pineapple/LR_stim/
addpath('~/OneDrive - University College London/1_DTMREIT/MRCI_toolbox')
VoxelSize = [3.55, 2.0833, 2.0833];
CurrAmpl = 0.001;       % injected current amplitude
Mask=int32(niftiread('Pinapple_mask_MREIT_RAS.nii.gz'));
Bz=niftiread('LR_meanBz_RAS.nii.gz');
J0=niftiread('LR_J0image.nii');
% J0(:,:,:,1)=(J0file(1).dat);
% J0(:,:,:,2)=(J0file(2).dat);
% J0(:,:,:,3)=(J0file(3).dat);
Bz(isnan(Bz))=0;
J0(isnan(J0))=0;
Mask(isnan(Mask))=0;
Mask=double(Mask);
Bz=double(Bz);
szJ0=size(J0);
szBz=size(Bz);
add=szBz(2)-szBz(1);
Bz=[Bz; zeros(add,szBz(2),szBz(3))];
J0=[J0; zeros(add,szJ0(2),szJ0(3), szJ0(4))];
Mask=[Mask; zeros(add,szBz(2),szBz(3))];



for thisSlice=1:szBz(3)
recon_parameters = reconstruction_parameters('projected current density type-1',...
    'VoxelSize',VoxelSize,'Mask',Mask(:,:,thisSlice));
thisJp=mrci_projected_current_density1(squeeze(J0(:,:,thisSlice,:)),squeeze(Bz(:,:,thisSlice)),recon_parameters) ;  
Jp(:,:,thisSlice,:)=thisJp;
disp(strcat('Just finished slice... ',num2str(thisSlice)))
end     

info=niftiinfo('LR_J0image.nii');
info.Filename=('/Users/ainsliej/OneDrive - University College London/1_DTMREIT/Pineapple/LR_stim/LR_Jp.nii.gz');
info.ImageSize=[84,84,96,3];
niftiwrite(Jp, 'LR_Jp', info);



%% now run the DTI bit 

clear all 
cd /Users/ainsliej/'OneDrive - University College London'/1_DTMREIT/Pineapple/
DTI=niftiread('DTI/Pinapple_tensor_resampled_RAS.nii.gz');
szDTI=size(DTI);
add=szDTI(2)-szDTI(1);
DTI=[DTI; zeros(add,szDTI(2),szDTI(3), szDTI(4))];
Mask=double(niftiread('LR_stim/Pinapple_mask_MREIT_RAS.nii.gz'));
Mask=[Mask; zeros(add,szDTI(2),szDTI(3))];
VoxelSize = [3.55, 2.0833, 2.0833];
CurrAmpl = 0.002;       % injected current amplitude

Jp_LR= niftiread('LR_stim/LR_Jp.nii');
Jp_AP= niftiread('AP_stim/AP2LR_Jp.nii.gz');
Mag_info=niftiinfo('AP_stim/AP2LR_Jp.nii.gz');
Mag_info.ImageSize=[84,84,96];
Mag_info.PixelDimensions=[3.55, 2.0833, 2.0833];

Jp_LR_mag=(Jp_LR(:,:,:,1).^2+Jp_LR(:,:,:,2).^2+Jp_LR(:,:,:,3).^2).^(1/2);
Jp_AP_mag=(Jp_AP(:,:,:,1).^2+Jp_AP(:,:,:,2).^2+Jp_AP(:,:,:,3).^2).^(1/2);
Mag_info.Filename=('/Users/ainsliej/OneDrive - University College London/1_DTMREIT/LR_stim/Jp_mag.nii');
niftiwrite(Jp_LR_mag, 'LR_stim/Jp_mag.nii', Mag_info);
Mag_info.Filename=('/Users/ainsliej/OneDrive - University College London/1_DTMREIT/AP_stim/Jp_mag.nii');
niftiwrite(Jp_AP_mag, 'AP_stim/Jp_mag.nii', Mag_info);

All_Jp=cat(5,Jp_LR, Jp_AP);
szJp=size(All_Jp);

for thisSlice=1:szJp(3)
recon_parameters = reconstruction_parameters('dtmreit noniterative',...
    'VoxelSize',VoxelSize,'Mask',Mask(:,:,thisSlice), 'Lambda', 0);
[thisC,scaleFactor]=mrci_dtmreit_noniterative(squeeze(All_Jp(:,:,thisSlice,:,:)),squeeze(DTI(:,:,thisSlice,:)),recon_parameters) ;  
Cond(:,:,thisSlice,:)=thisC;
ScaleF(:,:,thisSlice,:)=scaleFactor;
disp(strcat('Just finished slice... ',num2str(thisSlice)))
end   


info=niftiinfo('LR_stim/Jp_2dir.nii.gz');
info.Filename=('/Users/ainsliej/OneDrive - University College London/1_DTMREIT/Pineapple/conductance.nii.gz');
info.ImageSize=[84,84,96,6];
info.Datatype='double';
niftiwrite(Cond, 'conductance', info);
info.Filename=('/Users/ainsliej/OneDrive - University College London/1_DTMREIT/Pineapple/ScaleFactor.nii.gz');
info.ImageSize=[84,84,96];
info.PixelDimensions=[3.55, 2.0833, 2.0833];
info.Datatype='double';
niftiwrite(ScaleF, 'ScaleFactor', info);


%% Visualising things

D=Mask; % which image do I want to see 

figure
contourslice(D,[],[],[1:52])
view(3)

  
%% Fit curve

function T2star= curvefit(TE, allEchos) 
[xData, yData] = prepareCurveData( TE, allEchos );

% Set up fittype and options.
ft = fittype( 'S0*exp(-x/T2s)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0];
opts.StartPoint = [60 100];
opts.Upper = [Inf 150];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
coefs= coeffvalues(fitresult);
T2star=coefs(2);
end 


