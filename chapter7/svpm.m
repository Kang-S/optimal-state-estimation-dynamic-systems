function [svm,obsvm] = svpm(lambda,bm,cm,n_markov)
%
%	   Computes modal observability matrix and singular values of the
%    modal participation to the pulse response samples.
%                                                              
%     function [sv,obsvm] = svpm(lambda,bm,cm,p)
%                                                              
%   INPUT PARAMETERS
%     lambda: n x 1 vector containing system eigenvalues                                 
%     bm    : n x r modal input matrix             
%     cm    : m x n modal output matrix                               
%     p     : number of Markov parameters                              
%                                                              
%   OUTPUT PARAMETERS
% 
%		   svm   : n x 1 vector containing the normalized singular values                                
%     obsvm : mp x n modal observability matrix                         

%         J. N. Juang 11-29-91
%         NASA Spacecraft Dynamics Branch                           
%                                                                   

%
%  Compute the modal observability matrix
%
  n=length(lambda);
      c0=cm;obsvm=c0;
      if n_markov > 1;
        for j=2:n_markov; 
           c0=c0*diag(lambda);
           obsvm=[obsvm;c0];
        end; 
      end;
%
%  Compute the sigular values of the modal participation to the
%  pulse response samples
%    
   svm=zeros(n,1);
    for j=1:n;
%     svm(j)=sqrt([obsvm(:,j)'*obsvm(:,j)]*[bm(j,:)*bm(j,:)']);
      svm(j)=sqrt([bm(j,:)*bm(j,:)']*[cm(:,j)'*cm(:,j)])/abs(1-abs(lambda(j)));
    end;
   svm=svm ./max(svm);
