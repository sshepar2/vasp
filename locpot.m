%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                   %%%%%%%%%%%%%%
%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%
%%%%%%%%%%%   %                                            %%%%%   %%%%%%%%
%%%%%%%%%%%   %                                              %%%%%   %%%%%%
%%%%%%%%%%%%   %                                             %%%%%%%   %%%%
%%%%%%%%%%%%   %                                             %%%%%%%%   %%%
%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%    %%%%%%%%%%%%            %%%%%%%%%   %%%
%%%%%%%%%%%%%%%    %%%%%%%         %%%%%%%%%%              %%%%%%%%%   %%%%
%%%%%%%%%%%%%%%%%%   %            %%%%%%%                 %%%%%%%%%   %%%%%
%%%%%%%%%%%%%%%%%   %            %%%%                    %%%%%%%%%   %%%%%%
%%%%%%%%%%%%%%%%   %            %                     %%%%%%%%%%%   %%%%%%%
%%%%%%%%%%%%%%%   %                                %%%%%%%%%%%   %%%%%%%%%%
%%%%%%%%%%%%%%   %                              %%%%%%%%%%%   %%%%%%%%%%%%%
%%%%%%%%%%%%%   %             %                    %%%%   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   %             %%%%%                  %%%%%%   %%%%%%%%%%%%%%
%%%%%%%%%%%   %             %%%%%%%                 %%%%%%%%   %%%%%%%%%%%%
%%%%%%%%%%   %             %%%%%%%%%               %%%%%%%%%%   %%%%%%%%%%%
%%%%%%%%%   %             %%%%%%%%%%              %%%%%%%%%%%   %%%%%%%%%%%
%%%%%%%%   %             %%%%%%%%%               %%%%%%%%%%%   %%%%%%%%%%%%
%%%%%%%   %                                     %%%%%%%%%%%   %%%%%%%%%%%%%
%%%%%%   %                                     %%%%%%%%%%%   %%%%%%%%%%%%%%
%%%%%   %                                     %%%%%%%%%%%   %%%%%%%%%%%%%%%
%%%%   %                                   %%%%%%%%%%%%%   %%%%%%%%%%%%%%%%
%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%
%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%
%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%
%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                            %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuart Shepard %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Binghamton University %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2015 - 2020 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VASP - LOCAL POTENTIAL PLOTTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script reads in the 3D local potential data by reading in a file
% structured like the LOCPOT output file from the VASP program. To make
% VASP produce a LOCPOT file, use either the tag LVTOT = True or LVHAR =
% True in your INCAR file. LVTOT gives the entire local potential while 
% LVHAR leaves out the contribution of the exchange correlation potential.
% 
% LVTOT = True ---> LOCPOT (Vions + Vh + Vxc)
% LVHAR = True ---> LOCPOT (Vions + Vh)
%
% Where Vh is the Hartree potential.
%
% See: https://www.vasp.at/wiki/index.php/LVHAR for more info and links.
%
% Currently (as of 03-18-2020), this script stores the xy-averaged local
% potential and the atomic positions in cartesian coodinates and plots the
% xy-averaged local potential versus z. One could keep the full 3D data and
% plot the potential data however you wish, as well as the atoms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% How to use %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  1. Enter file path to LOCPOT structured file in the 'Getting Preliminary 
%     Data' section.
%
%  2. Enter fermi energy manually in the 'Plotting' section. Fermi energy
%     can be found by using the command 'grep fermi OUTCAR' in your 
%     terminal in the directory where VASP produced the LOCPOT file.
%
%  3. Change plotting preferences as preferred in the 'Plotting' section.
%
%  4. Run the script.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
%% Getting Preliminary Data %%

% enter file name below
fid = fopen('LOCPOT');
line = fgetl(fid);
line = fgetl(fid);
line = fscanf(fid,'%f %f %f',3);
a1 = line;
line = fscanf(fid,'%f %f %f',3);
a2 = line;
line = fscanf(fid,'%f %f %f',3);
a3 = line;

A = [a1 a2 a3];

fgetl(fid);
line = fgetl(fid);
atoms = regexp(line, ' ', 'split'); atoms = atoms(~cellfun('isempty',atoms)); asize = size(atoms);
line = fscanf(fid,'%i',asize(2));
num_atoms = line;
tot_atoms = sum(num_atoms);

% getting type of coordinates
fgetl(fid);
coord_type = fgetl(fid);
lett = coord_type(1);

% checking if selective dynamics line is used
if lett == 'S' or lett == 's'
    coord_type = fgetl(fid);
    lett = coord_type(1);
end

% getting atomic_coordinates if wish to plot atoms in 3D
atom_coords = zeros(tot_atoms,3);
for m=1:tot_atoms
    line = fscanf(fid,'%f %f %f',3);
    atom_coords(m,:) = line;
end

% converting to cartesian coordinates
if lett == 'D'
    atom_coords = (A*atom_coords')';
end


fgetl(fid);
fgetl(fid);
grid = fscanf(fid,'%i %i %i',3);

%% Getting Potential Data (Vh, Vh+Vxc...) %%

pots = fscanf(fid,'%f');

potmat = zeros(grid(1),grid(2),grid(3));

for c=1:grid(3)
    for b=1:grid(2)
        for a=1:grid(1)
            potmat(a,b,c) = pots(a + (b-1)*grid(1) + (c-1)*grid(2)*grid(1));
        end
        
    end
    
end

clear pots;

% calculating z points for plotting
z = 1:grid(3);
z = z - 1;
z = z.*sqrt(a3(1)^2 + a3(2)^2 + a3(3)^2)./grid(3);

one = potmat(:,:,1);
one = one(:);

% calculating avegares over xy-planes
xy_ave = zeros(1,grid(3));
for c=1:grid(3)
    one = potmat(:,:,c);
    one = one(:);
    xy_ave(1,c) = mean(one);
    
end

clear potmat;

%% Plotting %%

% enter fermi energy manually or just set to zero
fermi = 0;

plot(z,xy_ave-fermi,'b');

% setting axis plotting limits
     %xmin %xmax  %ymin               %ymax
axis([z(1) z(end) min(xy_ave)-fermi-1 max(xy_ave)-fermi+1]);

% setting general properties
set(gca,'Fontsize',30,'FontName','Times New Roman');

% setting x axis properties
%xticks([num1 num2..]); % xtick values
%xticklabels({'label1','label2',..}); % xtick labels
xlabel('z (Angstrom)');

% setting y axis properties
%yticks([-25 -20 -15 -10 -5 0 5 10 15]); %ytick values
%yticklabels({'-25','-10',}); %ytick labels
ylabel('Potential E - E_F (eV)'); %yaxis label


pbaspect([1 1 1]); % change plot shape proportions
box on; % include axis ticks on both sides of plot
 
hold;
