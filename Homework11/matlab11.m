%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Program 11 Supplement: Matlab function file to plot the
%% u versus (x,y) data from program11.cpp
%%
%% Inputs:
%% N, M              Grid parameters for solution
%% program11.out     Data file from program11.cpp; first
%%                   few lines of text must be removed
%%                   from file before using.
%%
%% Outputs:
%% prog11_Graph.pdf  PDF file of graph of u(x,y)
%% prog11_Cmap.pdf   PDF file of contour map of u(x,y)
%%
%%
%% Here's how to get started:
%%
%%  1) Copy matlab11.m (this file) into your working directory.
%%
%%  2) Launch the software MATLAB by typing "matlab" at the
%%     Linux prompt in your working directory.  (The math 
%%     dept machines have this software installed.)  This
%%     may take a few seconds -- you should get a MATLAB
%%     window with a prompt that looks like ">>".
%%
%%  3) Adjust the values of N,M below.  Remove the first
%%     few lines of text from program11.out.
%%
%%  4) To run this program, type "matlab11" at the prompt in 
%%     the MATLAB window.
%%
%%  5) After the program runs, the results will be displayed
%%     on the screen and saved in the files "prog11_Graph.pdf" 
%%     and "prog11_Cmap.pdf" in your working directory.
%%
%%  6) To re-run the program, just type "matlab11" at the
%%     MATLAB prompt again.  If you are using a new data file,
%%     remember to remove the first few lines of text.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clear MATLAB workspace
clear all ;  

%% Specify size of grid
N=29; M=29;


%% Load program output into an array and
%% extract list of x-values and y-values
%% and table of u-values.
xyuData = importdata('program11.out') ;
for i=1:N+2
  xvec(i) = xyuData(1+(i-1)*(M+2),1) ;
end
for j=1:M+2
  yvec(j) = xyuData(j,2) ;
end
for i=1:N+2
  for j=1:M+2
    umat(i,j) = xyuData(j+(i-1)*(M+2),3) ;
  end
end

%% Make surface plot (can adjust view by 
%% dragging with mouse)
figure(1) ; clf ;
surf(xvec,yvec,umat') ;
xlabel('x') ;
ylabel('y') ;
zlabel('z') ;
title('Graph of u(x,y)') ;
grid on ;
rotate3d on ;

%% Print surface plot
print -dpdf prog11_Graph.pdf ;

%% Make contour map 
NumContourLines = 20 ;
figure(2) ; clf ;
contourf(xvec,yvec,umat',NumContourLines) ;
colorbar('location','eastoutside') ;
colormap(jet(NumContourLines)) ;
xlabel('x') ;
ylabel('y') ;
title('Contour map of u(x,y).  Color bar indicates u-values.') ;
grid on ;

%% Print contour map
print -dpdf prog11_Cmap.pdf ;
