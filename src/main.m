%% Initialization
   clear
   clc
   Main_Menu
   
%% Main Menu
function Main_Menu
clc
clear
disp('╔════════════════════════════════╗')
disp('║ Choose an option:                                  ║')
disp(['║   ([' 8 '1]' 8 ') Hohmann Transfer.                            ║'])  
disp(['║   ([' 8 '2]' 8 ') Bi-Elliptic Transfer                         ║'])  
disp(['║   ([' 8 '3]' 8 ') One-Tangent Transfer.                        ║'])
disp(['║   ([' 8 'X/x]' 8 ') to terminate.                              ║'])
disp('╚════════════════════════════════╝')
fprintf('-------------------------------------\n')
UserInput = input('~$ ','S'); %the dollar is inpired by Kali Linux terminal

    if UserInput == "1"
                Hohmann_Core();

            elseif UserInput == "2"    
                BiElliptic_Core();

            elseif UserInput == "3" 
                OneTangent_Core();
                
            elseif UserInput == "X" || UserInput == "x"  
              % If the user input is X/x, this will go to function EndNa
                EndNa;
        else    
                Main_Menu;
    end
end

% Hohmann Transfer Function

% Terminate
function EndNa
    clear
    clc
    fprintf(2, 'User Terminated \n')
    return;
end

