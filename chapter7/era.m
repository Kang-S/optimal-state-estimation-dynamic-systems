function [a,b,c,d,xs,xt,mac,p,q] = era(y,m,r,n,nm,dt)
%                                                              
%	   Construct a state space model from pulse respone samples
%    using the ERA time domain technique.
%                                                              
% [A,B,C,D,xs,xt,mac] = era(y,m,r,n,nm,dt)     
%                                                              
% y   : [y_1 y_2 y_3... y_ni];                                 
%       y_i = response samples due to a unit pulse            
%             at the i-th input                               
% m   : no. of system ouputs.                              
% r   : no. of system inputs.                              
% n   : desired model order                                 
%       (n=0; User's interaction is required to determine
%        the system order by looking at the singular value plot.)   
% nm  : no. of block row shifts to form H(0) (Hankel Matrix)                          
% dt  : data sample time interval                          
%                                                              
% The identified model is                                  
%	     x(k+1) = Ax(k) + Bu(k)                                
%   	  y(k) = Cx(k) + Du(k)                                
% xs  : singular values of the correlation matrix.  	     
% xt  : 1st column = Dampings (%)                         
%       2nd column = Frequencies (Hz)                     
%       3nd column = Continuous-time complex eigenvalues                    
% mac : 1st column = Singular value of the pulse response samples 
%                    contributed by the identified individual modes 
%       2nd column = Modal Amplitude Coherence (MAC)                  

%         J. N. Juang 11-29-91
%         NASA Spacecraft Dynamics Branch                           
%                                                                   
format short e; format compact
   [nd,i]=size(y);
   nr=nd-nm;% All the data will be used.
   flagera=1; if n==0; flagera=0; n=200; end;   
   if nr<1;
     disp(['Errors in the eradc m-file occur.  The data length is not'])
     disp(['long enough.  Try to increase the data length or reduce the'])
     disp(['number of row shifts nm. See User guide for more information.'])
     return;
   end; 
%
% Generate the Hankel matrix H:
%
    disp(['ERA is used now.']);
    disp(['The Hankel matrix size for ERA is ' num2str(nm*m) ' by '...
           num2str((nr-1)*r) '.']);
    disp(['If the size is too big for your computer to handle,']);
    disp(['it is better to use function eradc instead to save']);
    disp(['computer memory and computational time.']);

    Markov=p2m(y,r);
    [d,H]=hankl(Markov,r,nm);clear Markov
%
% Singular value decomposition:
%
    [p,xs,q]=svd(H(:,1:(nr-1)*r));
    xs=diag(xs);
     disp(['Maximum Hankel singular value = ' sprintf('%e',max(xs))])
     disp(['Minimum Hankel singular value = ' sprintf('%e',min(xs))])
%
% Determine the order of the system
%
    if flagera==0; 
      semilogy([xs],'*');
      xlabel(' Number');ylabel('SV Magnitude');
      title(' Hankel Matrix Singular Values');pause;
    end;
 sigcont=sum(xs);   
 while n > 0; 
   if flagera==0;
      n=input('Desired Model Order (0=stop)=: ');
      if isempty(n)==1;break;end;if n==0;break;end;
      sigkpt=sprintf('%g',100*sum(xs(1:n))/sigcont);
      disp([' Model Describes ' sigkpt ' (%) of Test Data'])
   end;
   if flagera==1;
       for i=1:nm*m;
          if (sum(xs(1:i,1))/sigcont>.9999999999);n_index=i;break;end;
       end;
       if n>n_index; 
         disp(['The initial order is set to ' num2str(n) '.'])
         disp(['It is now set to ' num2str(n_index) '.']) 
         n=n_index;
       end;       
     end;
%
% Calculate a discrete model realization
%
    d1=sqrt(xs(1:n)); 
    d2=1.0 ./d1;
    a=diag(d2)*p(:,1:n)'*H(:,r+1:nr*r)*q(:,1:n)*diag(d2);
    b=diag(d1)*q(1:r,1:n)';
    c=p(1:m,1:n)*diag(d1);

%
%  Sort the eigenvalues and eigenvectors
%
    [g,xf]=eig(a);
    [lambda,index]=sort(diag(xf));clear xf;
    g=g(:,index);
%
%   Calculate the singular values of modal participation to
%   the pulse response samples and modal amplitude coherence
%
   mac=zeros(n,2);
   [mac(:,1),cq]=svpm(lambda,g\b,c*g,nm);
   rq=p(:,1:n)*diag(d1)*g;clear g;
   mac(:,2)=abs(diag(rq'*cq) ./sqrt(diag(rq'*rq) .*diag(cq'*cq)));
   xt=deg2hz(lambda,dt);
%   if flagera==0;
     disp('   Damping(%)    Freq(HZ)      Mode SV       MAC')
	     disp([xt(:,1:2) mac])
%   end;
%
   if flagera==1; n=0; end;
 end;
