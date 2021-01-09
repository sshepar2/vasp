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
% Stuart Shepard %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Binghamton University %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2015 - 2020 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VASP - WAVE FUNCTION PLOTTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script is to be used along with the fortran script, WaveTransPlot,
% created by:
% 
% R. M. Feenstra and M. Widom
% Department of Physics, Carnegie Mellon University, Pittsburgh, PA 15213
%
% Info and downloads found at:
%
% https://www.andrew.cmu.edu/user/feenstra/wavetrans/     (as of Feb.2020)
% 
% The WaveTranPlot executable converts the WAVECAR file (a binary file, 
% not readable as text) output by VASP to a .txt file. Because the contents
% can be large, usually one converts one specific wave function (Psi_nk) at
% a time.  The command can be run as,
% 
% >> ./WaveTransPlot -f WAVECAR -s s -k k -b b (-x coord -y coord)
% 
% Where one specifies the spin (-s 1 or 2), k-point (-k see OUTCAR or 
% vasprun.xml for numbering of the k-points), and band (-b see VASP output
% files for numbering). See the link above for more usage info.
%
% So as an example:
%
% >> ./WaveTransPlot -f WAVECAR -s 1 -k 50 -b 1
% 
% will output a file GCOEFF.txt
% 
% The file structure is as follows:
% 
% b1x b1y b1z
% b2x b2y b2z
% b3x b3y b3z
% spin# k-point# band#      ::just like the input command
% k1 k2 k3                  ::in fractional coordinates of rec lat vectors
% (eigenvalue, 0.0 Im part?) occupation
% ng1 ng2 ng3 (Re(c_G), Im(c_G))  ::Rec. lat. vec. indices, complex coeffs.
% ...
% ...
% repeats for spinor 2 (if spin-orbit coupling) :: No heading inbetween
% 
% The present Matlab script (psinks.m) is written to take in a file of this
% type (a single wave function). This script will not work for a
% WAVECAR output from a VASP spin-orbit-coupling calculation. Though this 
% script can be altered fairly straightforwardly (with the use of the 
% WaveTransSpinor executable, whose output .txt has the same structure as 
% GCOEFF.txt but with a second set of coefficients for the second spinor 
% following all the coefficients for the first spinor) to plot spinor 
% wavefunctions as well.
%
%
% The Kohn-Sham (KS) orbitals have the form:
%
% psi_nk(r) = 1/V^(1/2) * sum_over_G [c_G_nk*exp(i(G+k)r)] 
%
% where V is the real space volume, G is a reciprocal lattice vector where 
% the sum is cut off according to the ENCUT tag in the INCAR file. 
% c_G_nk is the complex planewave coefficient for the reciprocal lattice 
% vector G, for KS wavefunction with wavevector k and band index n.
%
%
%
% Using psink.m
%
% 1. Replace the path to the .txt you generated using WaveTransPlot
% in the 'Retrieving Data from generic GCOEFF.txt' section.
% 
% 2. Only if reading a .txt for a spinor (spin-orbit coupling,
% non-collinear spin), make the necessary changes in the 'Manipulating
% Data' section. (see comments)
%
% 3. Change isosurface values in the 'Plotting Wave Function' section, and
% edit plotting specs if necessary.
%
% Future Improvement List
% 1. Plot a real space unit cell (no more input files required)
% 2. Plot atom locations (require a POSCAR or some other VASP output file)
%
% Other Notes:
%
% VASP wave functions are pseudo-wave functions and do not represent the
% true real space wave function near the atomic nucleus.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Retrieving Data from generic GCOEFF.txt %%

clear;

% Replace with path to file.txt here
fid = fopen('GC1501.txt');
line = fscanf(fid,'%f %f %f',3);
b1 = line;
line = fscanf(fid,'%f %f %f',3);
b2 = line;
line = fscanf(fid,'%f %f %f',3);
b3 = line;

line = fscanf(fid,'%i %i %i',3);
skb = line; %[spin comp, k-point #, band #]

line = fscanf(fid,'%f %f %f',3);
k = line; % k-point in fractional coords of reciprocal lattice vectors

line = fgetl(fid);
line = fgetl(fid);
line = strsplit(line,{'(',',',')'});
line = cellfun(@str2num,line,'un',0);
eig = line{1,2}; % wave function eigenvalue
occ = line{1,4}; % occupation of wave function

% Retrieving reciprocal lattice vector indices and real and imaginary part
% of the planewave coefficients
coeffmat = textscan(fid,'%d%d%d%f%f','Delimiter',',()');
pw = size(coeffmat{1,1});
pw = pw(1); % total number of planewaves describing this KS orbital
% If the wave function is a spinor simply divide pw by 2 to get the total
% number of planewaves

%% Manipulating Data %%

% Expressing k-point in x,y,z components
kx = k(1)*b1(1) + k(2)*b2(1) + k(3)*b3(1);
ky = k(1)*b1(2) + k(2)*b2(2) + k(3)*b3(2);
kz = k(1)*b1(3) + k(2)*b2(3) + k(3)*b3(3);
kmag = sqrt(kx^2 + ky^2 + kz^2); % spatial frequency 2pi radians
lambda = 2*pi()/kmag; % wave function wave length
lx = 2*pi()/kx;
ly = 2*pi()/ky;
lz = 2*pi()/kz;
l3 = lx*ly*lz;

% Real space values at which the wave function will be calculated
% Can be whatever values you want. But! 
% Make x,y,z range at least larger than lx,ly,lz respectively so that you 
% can view at least one full k modulation of the wavefunction
x = (-3:3/15:20);
y = (-3:3/15:20);
z = (-8:3/15:20);
[X,Y,Z] = meshgrid(x,y,z);

% Expressing G + k in x,y,z components
Gk = zeros(pw,3);
for n=1:pw % for second spinor change for loop range to n=pw+1:2*pw (I think)
    Gplusk = [double(coeffmat{1,1}(n)) + k(1), double(coeffmat{1,2}(n)) + k(2), double(coeffmat{1,3}(n)) + k(3)];
    Gx = Gplusk(1)*b1(1) + Gplusk(2)*b2(1) + Gplusk(3)*b3(1);
    Gy = Gplusk(1)*b1(2) + Gplusk(2)*b2(2) + Gplusk(3)*b3(2);
    Gz = Gplusk(1)*b1(3) + Gplusk(2)*b2(3) + Gplusk(3)*b3(3);
    Gk(n,:) = [Gx, Gy, Gz];
end

% Building total wave function by summing over all periodic functions
% weighted by their coefficents.
psink = 0.*X;
for m=1:pw % for second spinor change for loop range to n=pw+1:2*pw (I think)
    psink = psink + (coeffmat{1,4}(m)+1i*coeffmat{1,5}(m)).*exp(1i.*(Gk(m,1).*X + Gk(m,2).*Y + Gk(m,3).*Z));

end

%% Plotting Wave Function %%

% Give system name for plot, optional.
SysName = 'BeS';

% Separating Real and Imaginary parts of the wave function for plotting
Repsink = real(psink);
Impsink = imag(psink);

% Plotting the wave function isosurface at equal and opposite maginitude
figure

% Negative Isosurface
isovalue1 = -0.2; % change to see different wave function isosurface
surf1 = isosurface(x,y,z,Repsink,isovalue1);
p1 = patch(surf1);
isonormals(x,y,z,Repsink,p1);
set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
daspect([1,1,1]);
view(3); axis tight
camlight; lighting gouraud

% Positive Isosurface
isovalue2 = 0.2; % change to see different wave function isosurface
surf2 = isosurface(x,y,z,Repsink,isovalue2);
p2 = patch(surf2);
isonormals(x,y,z,Repsink,p2);
set(p2,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface

hold;

% Plotting wave function wave vector direction and wavelength for guidance.
plot3([(max(x)+min(x))/2-lambda*kx/kmag/2,lambda*kx/kmag/2+(max(x)+min(x))/2],[(max(y)+min(y))/2-lambda*ky/kmag/2,lambda*ky/kmag/2+(max(y)+min(y))/2],[(max(z)+min(z))/2-lambda*kz/kmag/2,lambda*kz/kmag/2+(max(z)+min(z))/2],'k','Linewidth',3);

% Plotting Details
xlabel('x');
ylabel('y');
zlabel('z');

% Getting fraction values if k point for plotting
cutoff = 1E-7; % k component is zero if less than cutoff
k1 = k(1); k2 = k(2); k3 = k(3);
k1(k1<cutoff)=0; k2(k2<cutoff)=0; k3(k3<cutoff)=0;% setting to zero
if k1 ~= 0
    [num, dem] = rat(k1);
    kb1 = sprintf('%d/%d',num,dem);
else
    kb1 = '0';
end
if k2 ~= 0
    [num, dem] = rat(k2);
    kb2 = sprintf('%d/%d',num,dem);
else
    kb2 = '0';
end
if k3 ~= 0
    [num, dem] = rat(k3);
    kb3 = sprintf('%d/%d',num,dem);
else
    kb3 = '0';
end

t = sprintf('System: %s \n Wave Function, \x3c8_{nk\x3c3} \n n: %d \t k: #%d, %s, %s, %s \t \x3c3: %d\n E: %d eV \t Occ: %d \n Plot Info \n iso: +/- %d',SysName,skb(3),skb(2),kb1,kb2,kb3,skb(1),eig,occ,abs(isovalue1));
%t = strcat('BeS',';','n=',num2str(skb(3)),';','k=',num2str(skb(2)),';','spin=',num2str(skb(1)),';','E=',num2str(eig),';','occupation=',num2str(occ),';','iso=+/-',num2str(abs(isovalue1)));
title(t) 



