function [sys,x0,str,ts] = mpp(t,x,u,flag)
%Perturbar observar

switch flag
  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  % Initialize the states, sample times, and state ordering strings.
  case 0
    [sys,x0,str,ts]=mdlInitializeSizes;

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%    
  % Return the outputs of the S-function block.
  case 3
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  % There are no termination tasks (flag=9) to be handled.
  % Also, there are no continuous or discrete states,
  % so flags 1,2, and 4 are not used, so return an emptyu
  % matrix 
  case { 1, 2, 4, 9 }
    sys=[];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Unexpected flags (error handling)%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Return an error message for unhandled flag values.
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end

% end timestwo

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts] = mdlInitializeSizes()

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;  %%%DEFINIR AQUI EL NUMERO DE SALIDAS DEL BLOQUE
sizes.NumInputs      = 2;  %%%DEFINIR AQUI EL NUMERO DE ENTRADAS DEL BLOQUE
sizes.DirFeedthrough = 1;   % has direct feedthrough
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
str = [];
x0  = [];
ts  = [250e-3 0];   % inherited sample time
%global tm;
%ts  = [tm 0];   % inherited sample time

% end mdlInitializeSizes

%
%=============================================================================
% mdlOutputs
% Return the output vector for the S-function
%=============================================================================
%
function sys = mdlOutputs(t,x,u)
%VARIABLES GLOVALES
global Vpv_anterior;
global Ipv_anterior;
global Ppv_anterior;
global Vref;
%Aqui se muestrean los valores de las señales de entrada
Vpv=u(1);
Ipv=u(2);
%entrada1=u(2);
Ppv= Vpv*Ipv;
 %ALGORITMO MPP
if (Ppv-Ppv_anterior)~=0
    if (Ppv-Ppv_anterior) >0
        if Vpv-Vpv_anterior >0
            Vref=Vref+1;
        else 
            Vref=Vref-1;
        end
    else
        if Vpv-Vpv_anterior >0
            Vref=Vref-1;
        else 
            Vref=Vref+1;
        
        end
    end
end


%REFRESH VARIABLES
Ppv_anterior=Ppv;
Vpv_anterior=Vpv;
Ipv_anterior=Ipv;
   



sal(1)=Vref;

sys = sal(1);

% end mdlOutputs

