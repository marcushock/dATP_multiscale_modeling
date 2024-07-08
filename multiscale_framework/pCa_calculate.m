%%%% Calculates Hill coefficient, ec50

% Adapted by Abby Teitgen based on code from: 

% Tewari, S. G., Bugenhagen, S. M., Palmer, B. M, Beard,
% D. A. (2016). Dynamics of corss-bridge cycling, ATP hydrolysis, force
% generation, and deformation in cardiac muscle. J. Mol. Cell Cardiol. 96:
% 11-25.

% Tewari, S. G., Bugenhagen, S. M., Vinnakota, K. C., Rice, J. J., Janssen,
% M. L., Beard, D. A. (2016). Influence of metabolic dysfuntion on cardiac
% mechanics in decompensated hypertrophy and heart failure. J. Mol.
% Cell Cardiol. 94: 162-175.

% Lopez, R. Marzban, B., Gao, X. Lauinger, E., Van den Bergh, F.,
% Whitesall, S. E., Converso-Baran, K., Burant, C. F., Michele, D. E.,
% Beard, D. A. (2020). Impaired myocardial energetics causes mechanical dysfunction
% in decompensated failing hearts. Function, 1(2): zqaa018

% Marzban, B., Lopez, R., Beard, D. A. (2020). Computational modeling of coupled
% energetics and mechanics in the rat ventricular myocardium. Physiome.


% Inputs: Two vectors of the same length, the first containing the concentration and
% the second the force. If doses of 0 are contained in the data,
% these are used as control values and the data are normalised by
% their mean value.

function [hillCoeff ec50]=pCa_calculate(X_pCa,Y_pCa)
% Deal with 0 dosage by using it to normalise the results
normalised=0;
if (sum(X_pCa(:)==0)>0)
    % Compute mean control response
    controlResponse=mean(Y_pCa(X_pCa==0));
    % Remove controls from dose/response curve
    Y_pCa=Y_pCa(X_pCa~=0)/controlResponse;
    X_pCa=X_pCa(X_pCa~=0);
    normalised=1;
end

% Hill equation sigmoid
sigmoid=@(beta,x)beta(2)+(beta(1)-beta(2))./(1+(x/beta(3)).^beta(4));
% cf_23(x) = D+(A-D)/(1+(x/C)^B)
% A is the lower asymptote so guess it with min(y)
% B is the Hill's slope so guess it with the slope of the line between first and last point.
% C is the inflection point (the concentration of analyte where you have
% half of the max response) so guess it finding the concentration whose
% response is nearest to the mid response.
% D is the upper asymptote so guess it with max(y)

% Calculate some rough guesses for initial parameters
minResponse=min(Y_pCa);
maxResponse=max(Y_pCa);
slope=(Y_pCa(end)-Y_pCa(1))/(X_pCa(end)-X_pCa(1));
[~,Idx]=min(abs((Y_pCa-((max(Y_pCa)-min(Y_pCa))/2))));

% Fit the curve and compute the values
[coeffs,r,J]=nlinfit(X_pCa,Y_pCa,sigmoid,[minResponse maxResponse X_pCa(Idx) sign(slope)]);
ec50=coeffs(3);
hillCoeff=coeffs(4);

