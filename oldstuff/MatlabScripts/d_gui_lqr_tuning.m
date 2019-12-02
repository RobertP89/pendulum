% D_GUI_LQR_TUNING
%
% D_GUI_LQR_TUNING creates an input dialog box for interactive LQR controller tuning. 
% In order to tune the LQR controller,
% the user is required to enter values for Q and R, the diagonal weighting matrices.
%
% Nomenclature:
% Q     Diagonal Matrix (semi-positive definite)
% R     Diagonal Matrix (positive definite)
%
% Copyright (C) 2003 Quanser Consulting Inc.
% Quanser Consulting Inc.

% input dialog box to set LQR controller gains
function [ Q, R ] = d_gui_lqr_tuning( )
global Gui_QR Gui_Trial_ID Gui_Cnt_Offset
% for internal use (needed at every iteration)
global fname Q_dim R_dim Offset_save prompt lines_nb resize Gui_QR_Prev

% initializations: only at the first-call
if ( Gui_Trial_ID == 0 )
    fname = [ 'lqr_tuning_logfile' '.txt' ];
    fid = fopen( fname, 'a' );
    fprintf( fid, '%% ############## START OF THE LQR TUNING TRIALS ############## \n' );
    fclose( fid );
    % known dimensions of the 'prompt' entities
    R_dim = 1; % a.k.a number of inputs (to set manually)
    save_dim = 1; % "wish to save?" : 1 element
    cnt_dim = 1; % "wish to continue?" : 1 element
    % final  offsets counts
    Q_dim = size( Gui_QR, 2 ) - R_dim - save_dim - cnt_dim;
    Offset_save = Q_dim + R_dim + 1;
    Gui_Cnt_Offset = Offset_save + save_dim;
    % define the GUI layout
    for i = 1 : Q_dim
        prompt{ i } = [ 'Enter Q(' num2str( i ) ',' num2str( i ) ') :'];
    end
    for i = 1 : R_dim
        prompt{ Q_dim + i } = [ 'Enter R(' num2str( i ) ',' num2str( i ) ') :'];
    end
    prompt{ Offset_save } = 'Save parameters to a text file? ( ''yes'' / ''no'' ) :';
    prompt{ Gui_Cnt_Offset } = 'Do you wish to continue? ( ''yes'' / ''no'' ) :';
    lines_nb = 1;
    resize = 'on';
    Gui_QR_Prev = Gui_QR;
end 

if strcmp( Gui_QR{ Gui_Cnt_Offset }, 'yes' )
    Gui_Trial_ID = Gui_Trial_ID + 1;
end

% calls the input dialog box
title = [ 'LQR Tuning Window' ' : Trial # ' num2str( Gui_Trial_ID ) ];
Gui_QR = inputdlg( prompt, title, lines_nb, Gui_QR, resize );

% process the user-input data (from the GUI)
if ~isempty( Gui_QR )
    % the "ok" button has been pressed and a non-empty cell array has been returned
    for i = 1 : Q_dim
        Q( i ) = str2num( Gui_QR{ i } );
    end
    for i = 1 : R_dim
        R( i ) = str2num( Gui_QR{ Q_dim + i } );
    end
    if strcmp( Gui_QR{ Offset_save }, 'yes' )
        % save the current trial's Q and R matrices to a text file (for future reference)
        % append the new information to the existing file
        fid = fopen( fname, 'a' );
        fprintf( fid, '%% Trial # %d :\n', Gui_Trial_ID );
        % first format:
        for i = 1 : Q_dim
            fprintf( fid, 'Q(%d,%d) = %.6G; ', i, i, Q( i ) );
        end
        fprintf( fid, '\n' );
        for i = 1 : R_dim
            fprintf( fid, 'R(%d,%d) = %.6G; ', i, i, R( i ) );
        end
        fprintf( fid, '\n' );
        % second format: the result can be copied and pasted to the design file
        % cell array format: 
        % Gui_QR_Prev = { 'Q(i,i)', [..], 'R(i,i)', [..], 'save_to_file', 'continue_gui_tuning' }
        fprintf( fid, 'Gui_QR = {' );
        for i = 1 : Q_dim
            fprintf( fid, ' ''%.6G'',', Q( i ) );
        end
        for i = 1 : R_dim
            fprintf( fid, ' ''%.6G'',', R( i ) );
        end
        fprintf( fid, ' ''no'', ''yes'' }; \n\n' );
        fprintf( '\nTrial # %d : Q & R have been saved to the file "%s".\n', Gui_Trial_ID , fname );
        fclose( fid );
        Gui_QR{ Offset_save } = 'no';
    end
    Gui_QR_Prev = Gui_QR;
else
    % "cancel" has been pressed
    % for consistency, re-initialize to the last, i.e. previously, 
    % validated, saved, and returned values
    Gui_QR = Gui_QR_Prev;
    for i = 1 : Q_dim
        Q( i ) = str2num( Gui_QR_Prev{ i } );
    end
    for i = 1 : R_dim
        R( i ) = str2num( Gui_QR_Prev{ Q_dim + i } );
    end
    % do not save anything to file
    Gui_QR{ Offset_save } = 'no';
    % also abort the iterative tuning loop
    Gui_QR{ Gui_Cnt_Offset } = 'no';
end
% convert vectors to diagonal matrices
Q = diag( Q );
R = diag( R );
% end of function 'd_gui_lqr_tuning( )'
