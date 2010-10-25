% moment generating function for 3 step initiation regulated model

function QI =  MGRI(lambda,kab,kba,k12,k21,k23)
QI = kab.*k12.*k23./(kba.^2.*lambda+kba.*k12.*lambda+kba.*lambda.*kab+2.*kba.*lambda.^2+k21.*kba.*lambda+k21.*lambda.*kab+k21.*lambda.^2+k23.*kba.*lambda+k23.*k12.*kab+k23.*k12.*lambda+k23.*lambda.*kab+k23.*lambda.^2+lambda.*k12.*kab+k12.*lambda.^2+lambda.^2.*kab+lambda.^3);

