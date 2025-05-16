%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB CODE FOR FIGURE 3 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('NMR_E_FOUR_MATS.mat');
%%% FIGURE 3 A raw NMR data %%%
figure;imagesc(NMR_E_RAW);colormap jet;colorbar; %caxis([-1, 1]);
%%% FIGURE 3 B reordered NMR data %%%
figure;imagesc(NMR_E_ALL);colormap jet;colorbar; %caxis([-1, 1]); 
%%% FIGURE 3 C sample version of interconnected communities %%%
figure;imagesc(NMR_E_170);colormap jet;colorbar; %caxis([-1, 1]); 
%%% FIGURE 3 D population version of interconnected communities %%%
figure;imagesc(NMR_E_170_5_EST);colormap jet;colorbar; %caxis([-1, 1]); 


