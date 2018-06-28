function con_eig=deg2hz(dis_eig,dt)
%
%     Compute dampings (%) and frequencies (Hz) from the       
%     eigenvalues of an identified discrete model.             
%                                                              
%        function con_eig=deg2hz(dis_eig,dt)                   
%                                                              
%     dis_eig:eigenvalue vector from an                        
%             identified discrete model.                       
%     dt     :sampling rate (sec)                              
%                                                              
%     con_eig:1st column = Dampings (%)                        
%             2nd column = Frequencies (Hz)                    
%             3nd column = continuous-time complex eigenvalues 
%

%         L.G. Horta  4-25-91                                       
%         NASA Spacecraft Dynamics Branch                           
%                                                                   
          lamda=log(dis_eig)/dt;
           frequency=abs(lamda)/(2*pi);
           damping=-100*real(lamda)./abs(lamda);
           con_eig=[damping frequency lamda];
%           disp('   Damp(%)  Freq(Hz)')
%           disp(con_eig)
