%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB CODE FOR FIGURE 1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FIGURE 1 A genomics data %%%
load('Spellman_FOUR_MATS.mat');
figure;imagesc(Spell_RAW);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 A Raw
figure;imagesc(Spell_ALL);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 A Reordered

%%% FIGURE 1 B proteomics data %%%
load('Proteomics_FOUR_MATS.mat');
figure;imagesc(PD_RAW);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 B Raw
figure;imagesc(PD_ALL);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 B Reordered

%%% FIGURE 1 C gene expression data %%%
load('TCGA_Kidney_Gene_R_ALL_2448_1777.mat');
figure;imagesc(TCGA_Kidney_Gene_R_RAW_2448_1777);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 C Raw
figure;imagesc(TCGA_Kidney_Gene_R_ALL_2448_1777);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 C Reordered

%%% FIGURE 1 D multi-omics data %%%
load('Seed_FOUR_MATS.mat');
figure;imagesc(Seed_RAW);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 D Raw
figure;imagesc(Seed_ALL);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 D Reordered

%%% FIGURE 1 E brain imaging data %%%
load('EPSI_FOUR_MATS.mat');
figure;imagesc(EPSI_RAW);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 E Raw
figure;imagesc(EPSI_ALL);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 E Reordered

%%% FIGURE 1 F environmental plasma metabolomics data %%%
load('Metabolites_FOUR_MATS.mat');
figure;imagesc(Metabolites_RAW);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 F Raw
figure;imagesc(Metabolites_ALL);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 1 F Reordered

%%% FIGURE 3 A B plasma metabolomics data %%%
load('NMR_E_FOUR_MATS.mat');
figure;imagesc(NMR_E_RAW);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 3 A Raw
figure;imagesc(NMR_E_ALL);colormap jet;colorbar; caxis([-1, 1]); %%% in Figure 3 B Reordered


