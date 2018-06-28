function Markov=p2m(pulse,ni)
%
%     Routine to Rearrange the pulse response in Markov        
%     parameters                                               
%                                                              
%        function Markov=p2m(pulse,ni)                         
%                                                              
%     pulse=[y_1 y_2 y_3...y_ni];                              
%            y_i = response samples due to a unit pulse        
%                  at the i-th input                           
%     Markov=[D CB CAB ... CA**PB]                             
%     ni      = Number of inputs                               

%         L. G. Horta and J. N. Juang  6-25-91
%         NASA Spacecraft Dynamics Branch                           

      [nrow,ncol]=size(pulse);
       no=fix(ncol/ni);  
     Markov=zeros(no,ni*nrow);
        for i=1:ni
           nskl=i:ni:ni*nrow;
           Markov(1:no,nskl)=pulse(:,(i-1)*no+1:i*no)';
        end

 
