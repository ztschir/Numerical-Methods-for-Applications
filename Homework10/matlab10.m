%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Program 10 Supplement: Matlab function file to plot the
%% temperature versus position data from program10.cpp
%%
%% Inputs:
%% program10.out     Data file from program10.cpp; first
%%                   few lines of text must be removed
%%                   from file before using.
%%
%% Outputs:
%% prog10fig.pdf     PDF file of temperature vs position
%%
%%
%% Here's how to get started:
%%
%%  1) Copy matlab10.m (this file) into your working directory.
%%
%%  2) Launch the software MATLAB by typing "matlab" at the
%%     Linux prompt in your working directory.  (The math 
%%     dept machines have this software installed.)  This
%%     may take a few seconds -- you should get a MATLAB
%%     window with a prompt that looks like ">>".
%%
%%  3) To run this program, type "matlab10" at the prompt in 
%%     the MATLAB window.
%%
%%  4) After the program runs, the results will be displayed
%%     on the screen and saved in the file "prog10fig.pdf" in 
%%     your working directory.
%%
%%  5) To re-run the program, just type "matlab10" at the
%%     MATLAB prompt again.  If you are using a new data file,
%%     remember to remove the first few lines of text.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clear MATLAB workspace
clear all ;  

%% Load data into an array
xyData = importdata('program10.out') ;

%% Make plot
figure(1) ; clf ;
xmin = -2 ; xmax = 2 ; ymin = 0 ; ymax = 600 ; 
plot(xyData(:,1),xyData(:,2),'r-','LineWidth',3) ;
axis([xmin xmax ymin ymax]) ;
set(gca,'XTick',[xmin:0.5:xmax]) ;
set(gca,'YTick',[ymin:50:ymax]) ;
xlabel('position') ;
ylabel('relative temperature') ;
grid on ;

%% Print plot
print -dpdf prog10fig.pdf ;

