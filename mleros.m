function varargout=mleros(Hx,Gx,thini,params,algo,bounds,aguess)
% [thhat,covh,lpars,scl,thini,params,Hk,k]=...
%          MLEROS(Hx,Gx,thini,params,algo,bounds,aguess)
%
% Performs a maximum-likelihood estimation for CORRELATED loads as in
% Olhede & Simons (2013) by minimization using FMINUNC/FMINCON.
%
% INPUT:
%
% Hx       Real matrix with surface and subsurface topography [m m],
% Gx       Real vector with Bouguer gravity anomaly data [m/s^2]
%          ... either from direct observation or from SIMULROS...
%          ... with the geographical coordinates linearly unwrapped...
% thini    An unscaled starting guess for the parameter vector with elements:
%          [D f2 r s2 nu rho] - see SIMULROS
% params   A parameter structure with constants assumed known, see SIMULROS
%          [DEL g z2 dydx NyNx blurs kiso] in the units of ...
%           kg/m^3 (2x), m/s^2, m (3x), "nothing" (3x), rad/m
%          blurs  0 Don't blur likelihood using the Fejer window
%                 N Blur likelihood using the [default: N=2] resampled Fejer window
%           kiso   wavenumber beyond which we are not considering the likelihood
% algo     'unc' uses FMINUNC for unconstrained optimization
%          'con' uses FMINCON with positivity constraints [default]
%          'klose' simply closes out a run that got stuck
% bounds    A cell array with those positivity constraints [defaulted]
% aguess    A parameter vector [s2 nu rho] that will be used in
%           simulations for demo purposes, and on which "thini" will be
%           based if that was left blank. If "aguess" is blank, there is
%           a default. If "thini" is set, there is no need for "aguess"
%
% OUTPUT:
%
% thhat    The maximum-likelihood estimate of the vector [scaled]:
%          [D f2 r s2 nu rho], in Nm, and "nothing", see SIMULROS
% covh     A Hessian-based covariance estimate of the parameters
% lpars    The logarithmic likelihood and its derivatives AT or NEAR the estimate
%          lpars{1} the numerical logarithmic likelihood [FMINUNC/FMINCON]
%          lpars{2} the numerical scaled gradient, or score [FMINUNC/FMINCON]
%          lpars{3} the numerical scaled second derivative, or Hessian [FMINUNC/FMINCON]
%          lpars{4} the exit flag of the FMINUNC/FMINCON procedure [bad if 0]
%          lpars{5} the output structure of the FMINUNC/FMINCON procedure
%          lpars{6} the options used by the FMINUNC/FMINCON procedure
%          lpars{7} any bounds used by the  FMINUNC/FMINCON procedure
%          lpars{8} the residual moment statistics used for model testing 
%          lpars{9} the predicted variance of lpars{8}(3) under the null hypothesis
% scl      The scaling applied as part of the optimization procedure
% thini    The scaled starting guess used in the optimization
% params   The known constants used inside, see above
% Hk       The spectral-domain interface topographies after deconvolution 
% k        The wavenumbers on which the estimate is actually based
%
% NOTE: 
%
% At least 'demo1' has been tested to run in an SPMD loop!
%
% EXAMPLE:
%
%% Perform a series of N simulations centered on th0
% mleros('demo1',N,th0)
%% Statistical study of a series of simulations done using 'demo1'
% mleros('demo2','02-Oct-2014')
%% Admittance/coherence study of a series of simulations
% mleros('demo3','02-Oct-2014')
%% Covariance study of a series of simulations
% mleros('demo4','02-Oct-2014')
%% One simulation and a chi-squared plot
% mleros('demo5')
% Attempting to find correlations where there is none
% mleros('demo6')
% Attempting to find the density contrast and compensation level
% mleros('demo7')
%
% Last modified by fjsimons-at-alum.mit.edu, 06/23/2026

% NOTE: There are demonstrably bad solutions when r is close to 0 and f
% is close to -1 and 1. Fix in bounding? Ignore, fix later? Also, keep
% the D from being unrealistically low. If rho is zero poops out.

if ~isstr(Hx)
  defval('algo','unc')
  if strcmp(algo,'con')
    % Parameters for FMINCON in case that's what's being used
    bounds={[],[],... % Linear inequalities
	    [],[],... % Linear equalities
	    [1e17 eps -0.99/100  eps 0.95/100  10]*100,... % Lower bounds
	    [Inf  10  +0.99      Inf 4.00     Inf],... % Upper bounds
	    []}; % Nonlinear (in)equalities
  else
    bounds=[];
  end

  % The necessary strings for formatting
  str0='%27s';
  str1='%12.0e ';
  str2='%12.5g ';

  % Supply the needed parameters, keep the givens, extract to variables
  fields={'DEL','g','z2','dydx','NyNx','blurs','kiso'};
  defstruct('params',fields,...
	    {[2670 630],9.81,35000,[20 20]*1e3,sqrt(length(Hx))*[1 1],2,NaN});
  struct2var(params)

  % The gravitational constant (in m^3/kg/s2)
  G=fralmanac('GravCst');

  % Being extra careful or not?
  defval('xver',0)
  
  % The parameters used in the simulation for demos, or upon which to base "thini"
  defval('aguess',[7e22 0.4 -0.75 var(Hx) 3.0 sqrt(prod(dydx.*NyNx))/5]);
  % Scale the parameters by this factor; fix it unless "thini" is supplied
  defval('scl',10.^round(log10(abs(aguess))));

  % Unless you supply an initial value, construct one from "aguess" by perturbation
  nperturb=0.25;
  % So not all the initialization points are the same!!
  defval('thini',abs((1+nperturb*randn(size(aguess))).*aguess))
  % Starting guess for correlation coefficient needs to be different
  if abs(thini(3))>1
    thini(3)=rand*2-1;
  end

  %erin add below, not sure if it is necessary yet, there is a problem with
  %scl being empty in demo2
  if ~isempty(inputname(2)) || any(aguess~=thini)
    scl=10.^round(log10(abs(thini)));
    disp(sprintf(sprintf('\n%s : %s ',str0,repmat(str1,size(scl))),...
		 'Scaling',scl))
  end

  disp(sprintf(sprintf('%s : %s ',str0,repmat(str2,size(thini))),...
	       'Starting theta',thini))
  
  % If you brought in your own initial guess, need an appropriate scale
  if ~isempty(inputname(3)) || any(aguess~=thini)
    scl=10.^round(log10(abs(thini)));
    disp(sprintf(sprintf('%s : %s ',str0,repmat(str1,size(scl))),...
		 'Scaling',scl))
  end
  % Now scale so the minimization doesn't get into trouble - bounds also
  thini=thini./scl;
    
  defval('taper',0)
  if taper==1
    % Were going to want to make a 2D taper - any taper
    disp(sprintf('%s with TAPERING, DO NOT DO THIS YET',upper(mfilename)))
    NW=2;
    E1=dpss(NyNx(1),NW,1);
    E2=dpss(NyNx(2),NW,1);
    Tx=repmat(E1,1,NyNx(2)).*repmat(E2',NyNx(1),1);
    % But should still watch out for the spectral gain I suppose, this isn't
    % done just yet, although it appears already properly normalized
    % However, this looks better now, doesn't it?
    Tx=Tx*sqrt(prod(NyNx));
    % Not doing anything still amounts to saying Tx=1
  else
    Tx=1;
  end

  % Create the appropriate wavenumber axis
  k=knums(params);
  
  % Modify to demean
  disp('NOT DEMEAN BOTH DATA SETS')
  Hx(:,1)=Hx(:,1)-mean(Hx(:,1));
  Gx=Gx-mean(Gx);
  % Let us NOT demean and see where we end up...

  % Turn the observation vector to the spectral domain
  Hk(:,1)=tospec(Tx(:).*Hx(:,1),params);
  Gk     =tospec(Tx(:).*Gx     ,params);
  Hk(:,2)=Gk.*exp(k(:).*z2)/2/pi/G/DEL(2);

  % Am I doing this right if I give it a new z2? Make sure. Keep the Hk
  % for output.
  
  if xver==1
    % This only if you're using the right z2! and no blurring
    % This is reliant on the gravity data to be amenable to deconvolution,
    % while we could fake it in the simulations by working with Hx. Check
    % quickly, and note that there are roundoff errors right away!
    % Is the normalization right? I recently absorbed this into TOSPEC.
    difer([tospec(Tx(:).*Hx(:,2),params)-Hk(:,2)]/length(Hk),8,[],NaN)
    % Should also compare this with what actually can come out of SIMULROS
    % itself, although with real data of course we don't have this.
  end
  
  NN=200;
  % And now get going with the likelihood using Hk(:,1:2) or [Hk(:,1) Gk]
  % [ off|iter|iter-detailed|notify|notify-detailed|final|final-detailed ] 
  % Should probably make the tolerances relative to the number of k points
  options=optimset('GradObj','off','Display','off',...
		   'TolFun',1e-11,'TolX',1e-11,'MaxIter',NN,...
		   'LargeScale','off');
  % The 'LargeScale' option goes straight for the line search when the
  % gradient is NOT being supplied.

  % Set the parallel option to (never) use it for the actual optimization
  % Doesn't seem to do much when we supply our own gradient
  options.UseParallel='always';

  if blurs==0 || blurs==1
    % Use the analytical gradient in the optimization, rarely a good idea
    % options.GradObj='on';
    if xver==1
      % Definitely good to check this once in a while
      options.DerivativeCheck='on';
    end
  end

  % And find the MLE! Work on scaled parameters
  try
      switch algo
        case 'unc'
          % disp('Using FMINUNC for unconstrained optimization of LOGLIROS')
          t0=clock;
          [thhat,logli,eflag,oput,grd,hes]=...
	      fminunc(@(theta) logliros(theta,params,Hk,k,scl),...
		      thini,options);
          ts=etime(clock,t0);
          % Could here compare to our own estimates of grad and hes!
        case 'con'
          % New for FMINCON
          options.Algorithm='active-set';
          % disp('Using FMINCON for constrained optimization of LOGLIROS')
          t0=clock;
          [thhat,logli,eflag,oput,lmd,grd,hes]=...
	      fmincon(@(theta) logliros(theta,params,Hk,k,scl),...
		      thini,...
      		      bounds{1},bounds{2},bounds{3},bounds{4},...
                      bounds{5}./scl,bounds{6}./scl,bounds{7},...
		      options);
          ts=etime(clock,t0);
        case 'klose'
          % Simply a "closing" run to return the options
          lpars{6}=options;
          lpars{7}=bounds;
          % Simply a "closing" run to return the options
          varargout=cellnan(nargout,1,1);
          varargout{end}=lpars;
          return
      end
      if xver==1
          disp(sprintf('%8.3gs per %i iterations or %8.3gs per %i function counts',...
		       ts/oput.iterations*100,100,ts/oput.funcCount*1000,1000))
      else
          disp(sprintf('\n'))
      end
  catch
    % If something went wrong, exit gracefully
    varargout=cellnan(nargout,1,1);
    return
  end
  
  % This is the entire-plane estimate (hence the factor 2!)
  covh=hes2cov(-hes,length(k(~~k))*2);

  % Talk!
  disp(sprintf(sprintf('\n%s : %s ',str0,repmat(str2,size(thhat))),...
	       'Estimated theta',thhat.*scl))
  disp(sprintf(sprintf('%s : %s\n ',str0,repmat(str2,size(thhat))),...
	       'Asymptotic stds',sqrt(diag(covh))))

  % Here we compute the moment parameters and recheck the likelihood
  [L, ~,momx]=logliros(scl.*thhat,params,Hk,k,scl);
  diferm(L,logli)
  %logliros (as opposed to logliosl) doesn't have the ability to make vr

  % Likelihood attributes
  lpars{1}=logli;
  lpars{2}=grd;
  lpars{3}=hes;
  lpars{4}=eflag;
  lpars{5}=oput;
  lpars{6}=options;
  lpars{7}=bounds;
  lpars{8}=momx;
  %lpars{9}=vr;

  % Generate output as needed
  varns={thhat,covh,lpars,scl,thini,params,Hk,k};
  varargout=varns(1:nargout);
elseif strcmp(Hx,'demo1')
  % If you run this again on the same date, we'll just add to THINI and
  % THHAT but you will start with a blank THZERO. See 'demo2'
  % How many simulations? the second argument after the demo id
  defval('Gx',[]);
  N=Gx; clear Gx
  defval('N',500)
  more off
  % What th-parameter set? The THIRD argument after the demo id
  defval('thini',[]);
  % If there is no preference, then that's OK, it gets taken care of
  th0=thini; clear thini
  % What fixed-parameter set? The FOURTH argument after the demo id
  defval('params',[]);

  % The number of parameters to solve for
  np=6;
    
  % Open files and return format strings
  [fids,fmts,fmti]=osopen(np);
 
  % Do it!
  good=0; 
  % Initialize the average Hessian
  avH=zeros(np,np);

  % Set N to zero to simply close THZERO out
  for index=1:N
    % Simulate data from the same lithosphere, watch the blurring
    [Hx,Gx,th0,p,k]=simulros(th0,params);

    % Check the dimensions of space and spectrum are right
    difer(length(Hx)-length(k(:)),[],[],NaN)

    % Form the maximum-likelihood estimate, pass on the params, use th0
    % as the basis for the perturbed initial values. Remember hes is scaled.
    t0=clock;
    [thhat,covh,lpars,scl,thini,p]=mleros(Hx,Gx,[],p,[],[],th0);
    ts=etime(clock,t0);

    % Initialize the THZRO file... note that the bounds may change
    % between simulations, and only one gets recorded here
    if ~any(isnan(thhat)) && index==1 % && labindex==1
      oswzerob(fids(1),th0,p,lpars,fmts)
    end

    % If a model was found, keep the results, if not, they're all NaNs
    % Ignore the fact that it may be at the maximum number of iterations
    % e=1

    % IF NUMBER OF FUNCTION ITERATIONS IS TOO LOW DEFINITELY BAD
    itmin=0;

    % A measure of first-order optimality (which in this
    % unconstrained case is the infinity norm of the gradient at the
    % solution)  
    % FJS to update what it means to be good - should be in function of
    % the data size as more precision will be needed to navigate things
    % with smaller variance! At any rate, you want this not too low.
    optmin=Inf;

    % Maybe just print it and decide later? No longer e>0 as a condition.
    % e in many times is 0 even though the solution was clearly found, in
    % other words, this condition IS a way of stopping with the solution
    % Remember that the correlation coefficient can be negative or zero!
    % The HS is not always real, might be all the way from the BLUROS?
    try % Because if there are NaNs or not estimate it won't work
      % Maybe I'm too restrictive in throwing these out? Maybe the
      % Hessian can be slightly imaginary and I could still find thhat
      if isreal([lpars{1} lpars{2}']) ...
	    && all(thhat([1 2 4 5 6])>0) ...
	    && all(~isnan(thhat)) ...
	    && lpars{5}.iterations > itmin ...
	    && lpars{5}.firstorderopt < optmin
	good=good+1;
	% Build the average of the Hessians for printout later
	avH=avH+lpars{3}; %MLEOSL has this (divided by) ./[scl(:)*scl(:)'];
	% Reapply the scaling before writing it out
	fprintf(fids(2),fmts{1},thhat.*scl);
	fprintf(fids(3),fmts{1},thini.*scl);

	% Print the optimization results and diagnostics to a file 
        oswdiag(fids(4),fmts,lpars,thhat,thini,scl,ts,var(Hx),covh)
      end
    end
  end
  % If there was any success at all, finalize the THZERO file
  % If for some reason this didn't end well, do an N==0 run.

  % Initialize if all you want is to close the file
  if N==0
    [Hx,Gx,th0,p,k]=simulros(th0,params); 
    good=1; avH=avH+1; 
    [~,~,lpars]=mleros(Hx,Gx,[],[],'klose');
    oswzerob(fids(1),th0,p,lpars,fmts)
  end
  
  if good>=1 
    % This is the average of the Hessians, should be close to the Fisher
    avH=avH/good;

    % Now compute the theoretical covariance and scaled Fisher
    % This we don't really need any other time, it's just for use to be
    % able to compare to the average of the Hessians in this simulation
    % PUT IN ABS FOR THE CORRELATION COEFFICIENT WHICH MAY BE NEGATIVE
    sclth0=10.^round(log10(abs(th0)));
    [covF,F]=covthros(th0./sclth0,p,k,sclth0);
    % Of course when we don't have the truth we'll build the covariance
    % from the single estimate that we have just obtained. This
    % covariance would then be the only thing we'd have to save.
    % if labindex==1
        oswzeroe(fids(1),sclth0,avH,good,F,covF,fmti)
    % end
  end

  % Put both of these also into the thzro file 
  fclose('all');
elseif strcmp(Hx,'demo2')
  defval('Gx',[]);
  datum=Gx;
  defval('datum',date)

  % The number of parameters to solve for
  np=6;

  % Load everything you know about this simulation
  [th0,thhats,params,truecov,covavhs,thpix,~,~,~,~,momx]=osload(datum); %mleosl also has [covXpix,covF0] at the end as output
  %[th0,thhats,params,truecov,E,v,~,~,momx]=osload(datum);

  % Report the findings of the moment parameters
  disp(sprintf('m(m(Xk)) %f m(v(Xk)) %f m(magic) %s v(magic) %f',...
	      mean(momx),var(momx(:,end))))

  % Plot it all - perhaps some selection on optis?
  [ah,ha]=mleplos(thhats,th0,truecov,E,v,params,sprintf('MLEROS-%s',datum));

  % Print the figure! Don't forget the degs.pl script
  figna=figdisp([],sprintf('%s_%s',Hx,datum),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna)); 
  system(sprintf('rm -f %s.eps',figna)); 
elseif strcmp(Hx,'demo3')
  defval('Gx',[]);
  datum=Gx;
  defval('datum',date)
  
  % The number of parameters to solve for
  np=6;

  % Load everything you know about this simulation
  [th0,thhats,params,truecov,E,v]=osload(datum);

  % Plot it all: one admittance/coherence curve for every estimate
  [ah,ha]=admiplos(thhats(randi(length(thhats),100,1),:),th0,truecov,E,v,params,[],length(thhats));
  
  % Make the plot
  figna=figdisp([],sprintf('%s_%s',Hx,datum),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna));
elseif strcmp(Hx,'demo4')
  defval('Gx',[]);
  datum=Gx;
  defval('datum',date)
  
  % The number of parameters to solve for
  np=6;

  % Load everything you know about this simulation
  [th0,thhats,params,truecov,E,v,obscov,sclcov]=osload(datum);

  % Make the plot
  ah=covplos(2,sclcov,obscov,truecov,params,thhats,th0,E,v,'ver');

  figna=figdisp([],sprintf('%s_%s',Hx,datum),[],1);
  system(sprintf('epstopdf %s.eps',figna)); 
elseif strcmp(Hx,'demo5')  
  % What th-parameter set? The SECOND argument after the demo id
  defval('Gx',[]);
  % If there is no preference, then that's OK, it gets taken care of
  th0=Gx; clear Gx
  % What fixed-parameter set? The THIRD argument after the demo id
  defval('thini',[]);
  params=thini; clear thini

  % Figure name
  figna=sprintf('%s_%s_%s',mfilename,Hx,date);

  % Simulate data, watch the blurring, verify COLCHECK inside
  xver=1;
  [Hx,Gx,th0,p,k,Hk]=simulros(th0,params,xver);
  
  % Initialize, take defaulted inside MLEROS for now
  thini=[];

  % Perform the optimization, whatever the quality of the result
  [thhat,~,logli,thini,scl,p,e,o,gr,hs]=mleros(Hx,Gx,thini,p);

  % Take a look at the unblurred gradient purely for fun, they should be
  % so small as to be immaterial
  grobs=-nanmean(gammakros(k,thhat.*scl,p,Hk))';
  
  % Take a look at the unblurred theoretical covariance at the estimate,
  % to compare to the observed blurred Hessian; in the other demos we
  % compare how well this works after averaging
  [covthat,F]=covthros(thhat,p,k,scl);

  % Collect the theoretical covariance for the truth for the title
  covth=covthros(th0./scl,p,k,scl);
 
  % Take a look at the scaled Fisher to compare with the scaled Hessian  
  F;
  hs;
  grobs;
 
  % Take a look at the scaled covariances
  predcov=covthat./[scl(:)*scl(:)'];
  % And compare to the inverse of the scaled Hessians in the full plane
  % Watch as this will change with the anticipated halfplane changes
  obscov=inv(hs)/length(k(:))*2;
  % These two should compare reasonably well in the unblurred case, I
  % would have thought - but of course it's ignoring stochastic
  % variability. If we can use the Hessian for the uncertainty estimation
  % we can do this for the cases where we can't come up with a
  % theoretical covariance, not even an unblurred one. Check std's.
  % Maybe should do this a hundred times?
  disp(sprintf('%s\n',repmat('-',1,97)))
  disp('predicted and observed scaled standard deviations and their ratio')
  disp(sprintf([repmat('%6.3f  ',1,length(obscov)) '\n'],...
	       [sqrt(diag(predcov))' ; sqrt(diag(obscov))' ; ...
		sqrt(diag(predcov))'./sqrt(diag(obscov))']'))
  disp(repmat('-',1,97))
  % Talk again!
  [str0,str2]=osdisp(th0,p);
  disp(sprintf(sprintf('%s : %s ',str0,repmat(str2,size(thhat))),...
	       'Estimated theta',thhat.*scl))
  disp(repmat('-',1,97))
  
  % Young's modulus 
  defval('E',1.4e11);
  % Poisson's ratio
  defval('v',0.25);

  % Quick plot, but see OLHEDESIMONS5
  clf
  ah=krijetem(subnum(2,3)); delete(ah(4:6)); ah=ah(1:3);
  % Maybe we should show different covariances than the predicted ones??

  % Time to rerun LOGLIROS one last time at the solution
  [L,~,momx]=logliros(thhat,p,Hk,k,scl);

  % Better feed this to the next code, now it's redone inside
  mlechiplos(2,Hk,thhat,scl,p,ah,0,th0,covth,E,v);
elseif strcmp(Hx,'demo6')
  % Distribution of the likelihood ratios for the correlation coefficient
  addpath('../olhede3')
  
  % Set up experiment
  params.NyNx=[64 64];
  params.blurs=2;
  xver=1;

  % No correlation!
  th0=[1e24 0.8 0.0025 2 2e4];
  nperturb=0.25;

  N=1000;
  logliesr=zeros(N,13);
  
  for index=1:N
    thini1=(1+nperturb*randn(size(th0))).*th0;
    % As if there might be a correlation, let's start from 0
    thini2=[thini1(1) thini1(2) 0 thini1(3:5)];
    
    % No correlation!
    [Hx,Gx,th0,p,k,Hk]=simulros0(th0,params,xver);
    
    % No correlation!
    [thhat1,~,logli1,thinisc1,scl1,p1,e1,o1,gr1,hs1]=mleros0(Hx,Gx,thini1,p,'con');
    % Possible correlation!
    [thhat2,~,logli2,thinisc2,scl2,p2,e2,o2,gr2,hs2]=mleros(Hx,Gx,thini2,p,'con');

    % Maybe should change the cellnan in the main script
    if isnan(thhat1)
      thhat1=nan(size(thini1));
    end
    if isnan(thhat2)
      thhat2=nan(size(thini2));
    end
    
    % Speak
    osdisp(thhat1.*scl1,p)
    osdisp(thhat2.*scl2,p)
    
    % Keep everything [1 7] are the likelihoods
    logliesr(index,:)=[logli1 thhat1.*scl1 logli2 thhat2.*scl2];    
  end

  % Contract to the not-a-not-a-numbers, then save
  logliesr=logliesr(~any(isnan(logliesr),2),:)
  % Do a little windsorization of the obvious outliers
  logliesr=logliesr(logliesr(:,1)>prctile(logliesr(:,1),1.1),:);
  logliesr=logliesr(logliesr(:,7)>prctile(logliesr(:,7),1.1),:);
  
  save logliesr logliesr params xver th0 Hx Gx p k Hk scl1 scl2
  
  % Young's modulus 
  defval('E',1.4e11);
  % Poisson's ratio
  defval('v',0.25);

  % Collect the theoretical covariance for the truth for the title and
  % the overlay. Work with the uncorrelated model
  covth1=covthros0(th0./scl1,p,k,scl1);
  sclcv1=covth1./[diag(sqrt(covth1))*diag(sqrt(covth1))'];
  obscv1=cov(logliesr(:,[2:6]));
  obscv1=obscv1./[diag(sqrt(obscv1))*diag(sqrt(obscv1))'];

  % Now work with what it would be under the correlated model 
  covth2=covthros([th0(1) th0(2) 0 th0(3:5)]./scl2,p,k,scl2);
  sclcv2=covth2./[diag(sqrt(covth2))*diag(sqrt(covth2))'];
  obscv2=cov(logliesr(:,[8:13]));
  obscv2=obscv2./[diag(sqrt(obscv2))*diag(sqrt(obscv2))'];

  % Likelihood ratio test
  X=length(k(:))*(logliesr(:,1)-logliesr(:,7));
  % Only one parameter in the nested models is different
  df=1;
  
  % Make a bunch of figures illustrating what we just did
  % First the results of the estimation itself
  figure(1)
  mleplos(logliesr(:,[2:6]),th0,covth1,E,v,p,'demo6 - uncorrelated')  
  figure(2)
  mleplos(logliesr(:,[8:13]),[th0(1) th0(2) 0 th0(3:5)],covth2,E,v,p,'demo6 - correlated') 

  % Now look at the covariance plots
  figure(3)
  covplos([],sclcv1,obscv1,covth1,p,logliesr(:,[2:6]),th0,E,v,[])
  figure(4)
  covplos([],sclcv2,obscv2,covth2,p,logliesr(:,[8:13]),[th0(1) th0(2) 0 th0(3:5)],E,v,[])
  
  % Now look at the distribution of the loglihood difference
  figure(5)
  ah=krijetem(subnum(2,2)); delete(ah(3:4)); ah=ah(1:2);

  xstr='loglihood difference X (A79)';
  xll=[0 3*2*df];
  xlls=[xll(1):df:xll(2)];

  % The below is more or less copied from MLECHIPLOS, should maybe spin off
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(1))
  [bdens,c]=hist(X,5*round(log(length(X))));
  bdens=bdens/indeks(diff(c),1)/length(X);
  bb=bar(c,bdens,1);
  set(bb,'FaceC',grey)
  hold on
  refs=linspace(0,max(X),100);
  plot(refs,chi2pdf(refs,df),'Linew',1,'Color','k')
  hold off
  xlim(xll)
  xl(1)=xlabel(xstr); 
  yl(1)=ylabel('probability density');
  axis square

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(2))
  % Note that SOME people use a different parameterization (b vs 1/b)
  h=qqplot(X,ProbDistUnivParam('gamma',[df/2 2])); 
  axis square; box on
  set(h(1),'MarkerE','k')  
  set(h(3),'LineS','-','Color',grey)
  % Extend the line to the full axis
  hold on
  xh=get(h(3),'xdata');
  yh=get(h(3),'ydata');
  h(4)=plot([xh(2) xll(2)],...
	    [yh(2) yh(2)+[yh(2)-yh(1)]/[xh(2)-xh(1)]*[xll(2)-xh(2)]]);
  set(h(4),'LineS','-','Color',grey)
  hold off
  top(h(3),ah(2))
  delete(get(ah(2),'ylabel'));
  delete(get(ah(2),'title'));
  delete(get(ah(2),'xlabel'));
  xlim(xll); ylim(xll)
  xl(2)=xlabel('observed X');
  yl(2)=ylabel('predicted X');

  % Cosmetics
  fig2print(gcf,'landscape')
  longticks([ah])
  %set([cat(1,findobj('FontSize',10)); yl(:); xl(:)],'FontSize',12)
  movev([ah],-.1)
  serre(ah,[],'across')
  serre(ah,[],'across')

  t=ostitle(ah,p,[],length(logliesr(:,[2:6])));
  movev(t,.25)
  
  for index=1:5
    figure(index)
    figna=figdisp([],sprintf('%s_%i','demo6',index),[],1);
    system(sprintf('degs %s.eps',figna));
    system(sprintf('epstopdf %s.eps',figna)); 
    system(sprintf('rm -f %s.eps',figna)); 
  end
elseif  strcmp(Hx,'demo7')
  % Try to find the compensation level from the likelihood
  addpath('../olhede3')
  
  % Set up experiment
  params.NyNx=[64 64];
  % Do both blurred and unblurred experiments
  params.blurs=0;
  xver=1;
  
  % No correlation!
  th0=[1e24 0.8 0.0025 2 2e4];
  nperturb=0.25;
  nperturb=0.025

  N=51;
  logliesrt=zeros(N,13);
  
  z2=linspace(10000,60000,N);
  
  % The RIGHT one here is number 26
  for index=1:6:N
    thini1=(1+nperturb*randn(size(th0))).*th0;
    % As if there might be a correlation, let's start from 0
    % CANNOT START FROM 0 IN THE UNBLURRED CASE!
    thini2=[thini1(1) thini1(2) 0 thini1(3:5)];
    
    % No correlation in the simulation
    [Hx,Gx,th0,p,k,Hk]=simulros0(th0,params,xver);
    
    % Save this junk so inside I can take a look at it
    save keepfornow Hx Gx th0 p k  Hk

    % Now assign the compensation depth and see where this leads
    p.z2=z2(index);
    disp(sprintf('Trying compensation depth z2 = %i',round(z2(index))))
        
    % No correlation in the inversion!
    [thhat1,~,logli1,thinisc1,scl1,p1,e1,o1,gr1,hs1,Hk1]=mleros0(Hx,Gx,thini1,p,'con',-1);

    %pause
    
    % Possible correlation in the inversion!
    [thhat2,~,logli2,thinisc2,scl2,p2,e2,o2,gr2,hs2,Hk2]=mleros(Hx,Gx,thini2,p,'con');

    % pause
    
    % Take a look at the actual Hk used internally, we shouldn't be using
    % the one that we made with the correct compensation depth!
    
    % Should look at the chi-2 here at every turn and reject it somehow    

    % Now should look at the chi-squared to see if we can discriminate that
    clf
    ah=krijetem(subnum(2,3)); 
    % Maybe should change the cellnan in the main script
    if isnan(thhat1)
      thhat1=nan(size(thini1));
      delete(ah(1:3))
    else
      % Speak - but should update the strings if we don't mean "TRUE" but "estimated"
      osdisp(thhat1.*scl1,p)
      mlechiplos(3,Hk1,thhat1,scl1,p1,ah(1:3),0,th0,[],[],[]);
    end
    if isnan(thhat2)
      thhat2=nan(size(thini2));
      delete(ah(4:6))
    else
      % Speak - but should update the strings if we don't mean "TRUE" but "estimated"
      osdisp(thhat2.*scl2,p)
      % Now should look at the chi-squared to see if we can discriminate that
      [cb,~,t]=mlechiplos(2,Hk2,thhat2,scl2,p2,ah(4:6),0,[th0(1) th0(2) 0 th0(3:5)],[],[],[]);
      if isnan(thhat1)
	try
	  delete(ah(1:3))
	end
	  movev(t,1)
      else
	delete(t)
      end
      movev([ah(4:6) cb],.1)
      % Last-minute novelty
      movev([ah(4:6) cb],-.05)
    end
    
    % Print figure - if nonempty
    figna=figdisp([],sprintf('%s_%i','demo7',z2(index)),[],1);
    system(sprintf('degs %s.eps',figna));
    system(sprintf('epstopdf %s.eps',figna)); 
    system(sprintf('rm -f %s.eps',figna)); 
    
    % Keep everything [1 7] are the likelihoods
    logliesrt(index,:)=[logli1 thhat1.*scl1 logli2 thhat2.*scl2];    
    
    drawnow
  end
  
  % Add the depth back in
  logliesrt=[z2(:) logliesrt];

  disp('This demo is NOT ready for no-interactive mode!')

  % Quick save to not lose it
  save logliesrt logliesrt params xver th0 Hx Gx p k Hk scl1 scl2

  % Do a little windsorization of the obvious outliers
  logliesrt=logliesrt(logliesrt(:,2)>0,:);
  logliesrt=logliesrt(logliesrt(:,8)>0,:);

  % Perhaps work from all saved data
  % load logliesrt_demo7
  
  % Quick plot
  figure(1)
  clf
  plot(logliesrt(:,1)/1000,logliesrt(:,2),'ro')
  hold on
  plot(logliesrt(:,1)/1000,logliesrt(:,8),'bv')
  hold off
  grid on
  xlabel('depth to z2 [km]')
  ylabel('likelihoods')
  set(gca,'xtick',[10 35 60],'xlim',[10 60],'box','on')
  longticks(gca,2)
  fig2print(gcf,'portrait')
  
  % Print to file
  figna=figdisp([],sprintf('%s','demo7_1'),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna)); 
  system(sprintf('rm -f %s.eps',figna)); 

  % Looks like this won't be able to discriminate. Look at the
  % chi-squares to rule them out?
  
  % At any rate, do it again and save a second set.

  % Make some plots
  figure(2)
  xlabs={'D','f^2','r','\sigma^2','\nu','\rho'};
  clf
  for index=1:3
    ah(index)=subplot(2,3,index);
    % Don't do it if there is no correlation to plot
    if index~=3
      plot(logliesrt(:,1)/1e3,logliesrt(:,2+index)-th0(index),'ro')
    end
    hold on
    % Remember the true correlation is zero
    plot(logliesrt(:,1)/1e3,logliesrt(:,2+6+index)-th0(index)*(index~=3),'bv')
    hold off
    grid on
    xlabel('depth to z2 [km]')
    title(sprintf('%s-%s_0',xlabs{index},xlabs{index}))
    ylim(th0(index)*[-1.25 1.25])
    if index==3; ylim([-1 1]/10); end
        
    ah(index+3)=subplot(2,3,3+index);
    plot(logliesrt(:,1)/1e3,logliesrt(:,4+index)-th0(2+index),'ro')
    hold on
    plot(logliesrt(:,1)/1e3,logliesrt(:,4+7+index)-th0(2+index),'bv')
    hold off
    grid on
    xlabel('depth to z2 [km]')
    title(sprintf('%s-%s_0',xlabs{3+index},xlabs{3+index}))
    ylim(th0(2+index)*[-1.25 1.25])
    if index==3; ylim(th0(2+index)*[-0.5 0.5]); end
  end
  set(ah,'xtick',[10 35 60],'xlim',[10 60],'box','on')
  longticks(ah)
  fig2print(gcf,'landscape')
  
  % Print to file
  figna=figdisp([],sprintf('%s','demo7_2'),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna)); 
  system(sprintf('rm -f %s.eps',figna)); 
end 

