function varargout=Tros(k,th,params,phi,xi,pxm,pxp,xver)
% [invTo,detTo,Lo,To]=Tros(k,th,params,phi,xi,pxm,pxp,xver)
%
% Calculates functions of isotropic T for the CORRELATED initial-loading
% scenario with the primary spectra S11 etc being the initial-loading
% ones as in Olhede & Simons (2013)
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
% th       The parameter vector with THREE elements (rest ignored)
%          D   Isotropic flexural rigidity [Nm]
%          f2  The subsurface-to-surface initial loading ratio
%          r   The subsurface-to-surface initial loading correlation
% params   A structure with AT LEAST these constants that are known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
% phi      Optionally precalculated phi, see PHIOS
% xi       Optionally precalculated xi, see XIOS
% pxm      Optionally precalculated (phi*xi-1)
% pxp      Optionally precalculated (phi*xi+1)
%
% OUTPUT:
%
% invTo    A 3-column vector with all the wavenumbers unwrapped,
%          invTo={invTo[1,1](k) invTo[1,2](k) invTo[2,2](k)}
% detTo    The determinant of To, a column vector over the wavenumbers
% Lo       The Cholesky factorization of To, as the lower-left matrix 
%          Lo={Lo[1,1](k) Lo[2,1](k) Lo[2,2](k)}, where Lo[1,2]=0
% To       The actual matrix To, in the same format as invTo
%
% EXAMPLE:
%
% [~,~,th0,p,k]=simulros([],[],[],1);
% [invTo,detTo,Lo,To]=Tros(k,th0,p,[],[],[],[],1);
%
% Last modified by fjsimons-at-alum.mit.edu, 03/19/2014

defval('xver',0)

% Extract the parameters from the input
D=th(1);
f2=th(2);
r=th(3);

DEL=params.DEL;
g=params.g;

defval('phi',phios(k,D,DEL,g));
defval('xi',xios(k,D,DEL,g));
% Note that this has a zero at zero wavenumber
defval('pxm',(phi.*xi-1));
defval('pxp',pxm+2);

% Forcefully set f2 to a positive number even if it means a throw back
f2=abs(f2);
% Precompute some factors
f=sqrt(f2);
ddxi=[DEL(1)+DEL(2).*xi];

% The inverse of T; ignore warnings as Inf turns to NaN in HFORMOS 
warning off MATLAB:divideByZero
% Eq. (A7)
fax=DEL(1)^(-2)*ddxi.^2/f2./pxm.^2;
warning on MATLAB:divideByZero
% First the part without the correlation as in dTros
invT=[fax.*(                   1+f2*dpos(DEL,2,-2)*phi.^2) ...
      fax.*(dpos(DEL,-1,1)*xi   +f2*dpos(DEL,1,-1)*phi   ) ...
      fax.*(dpos(DEL,-2,2)*xi.^2+f2)];
% Then the part of the correlation, as in dTros
warning off MATLAB:divideByZero

% You still MIGHT have exactly zero... as in MLEROS('demo7')
if r~=0
  % Eq. (A9)
  fax=dpos(DEL,-1,-1)*ddxi.^2/r/f./pxm.^2;
  warning on MATLAB:divideByZero
  invDT=[fax.*(               2*phi) ...
	 fax.*(  dpos(DEL,-1 ,1)*pxp) ...
	 fax.*(2*dpos(DEL,-2,2)*xi )];
  % And now the complete thing
  % You STILL want an inverse even if there is no correlation, so you need
  % to force the issue - otherwise Matlab won't do 0*Inf
  invTo=[invT-r^2*invDT]/(1-r^2);
else
  % Don't ruin this by doing 0*Inf which turns it all to NaN
  invTo=invT;
end

if nargout>=2 || xver==1
  % Compute the determinant of To; this will be zero at k=0
  % Eq. (A6)
  detT= f2*DEL(1)^4./ddxi.^4.*pxm.^2;
  % Eq. (A8)
  detDT=-r^2*detT;
  % Eq. (A10)
  detTo=detT+detDT;
else
  detTo=NaN;
end

if nargout>=3 || xver==1
  % Compute the Cholesky factorization of T
  % Eq. (A4)
  fax=ddxi.^(-1)./...
      sqrt(DEL(2)^2*xi.^2+f2*DEL(1)^2-2*r*f*dpos(DEL,1,1).*xi);
  Lo=[fax.*(DEL(2)^2*xi.^2+f2*DEL(1)^2-2*r*f*dpos(DEL,1,1).*xi) ...
      fax.*(-dpos(DEL,1,-1).*[DEL(2)^2*xi+f2*DEL(1)^2.*phi]+r*f*DEL(1)^2.*pxp) ...
      fax.*f*DEL(1)^2.*pxm.*sqrt(1-r^2)];
else
  Lo=NaN;
end

if nargout>=4 || xver==1 
  % Compute T itself, which is required when producing blurred things
  % Eq. (A2)
  fax=DEL(2)^2./ddxi.^2;
  % First the part without the correlation
  T=[fax.*(             xi.^2+f2*dpos(DEL,2,-2)        ) ...
     fax.*(-dpos(DEL,1,-1)*xi-f2*dpos(DEL,3,-3)*phi    ) ...
     fax.*( dpos(DEL,2,-2)   +f2*dpos(DEL,4,-4)*phi.^2)];
  % Then the part of the correlation
  % Eq. (A3)
  fax=r*f*fax;
  DT=[fax.*(-2*dpos(DEL,1,-1).*xi ) ...
      fax.*(   dpos(DEL,2,-2).*pxp) ...
      fax.*(-2*dpos(DEL,3,-3).*phi)];
  % And now the complete thing
  % Eq. (A1)
  To=T+DT;
else
  To=NaN;
end

% And now for the output
varns={invTo,detTo,Lo,To};
varargout=varns(1:nargout);

% Verification mode
if xver==1
  disp('Tros being verified')
  % Explicit verification of the determinant
  detcheck(detT ,T ,8)
  detcheck(detDT,DT,8)
  detcheck(detTo,To,8)
  
  % Explicit verification of the inverse by the Cayley-Hamilton theorem
  invcheck(invT ,detT ,T ,3,1)
  invcheck(invDT,detDT,DT,3,1)
  invcheck(invTo,detTo,To,3,1)

  % Explicit verification of the inverse by checking the identity
  invcheck(invT, detT, T ,3,2)
  invcheck(invDT,detDT,DT,3,2)
  invcheck(invTo,detTo,To,3,2)
  
  % Check the Cholesky by multiplication
  cholcheck(Lo,To,5,1)

  % Check the Cholesky by factorization at a random wave number
  cholcheck(Lo,To,4,2)
end
