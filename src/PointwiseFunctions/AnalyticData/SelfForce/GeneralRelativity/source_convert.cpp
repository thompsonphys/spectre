int psiTommy(int m, double a, double r, double th, double realBarry[10],  double imagBarry[10], double realTommy[10], double imagTommy[10]){
double factor = r/(2*M_PI);
double rplus = 1 + sqrt(1 - a * a);
double rminus = 1 - sqrt(1 - a * a);
double mdphi = m * a / (rplus - rminus) * log((r - rplus) / (r - rminus));
double cosmdphi = cos(mdphi);
double sinmdphi = sin(mdphi);
realTommy[0] = factor*(realBarry[0]*cosmdphi+imagBarry[0]*sinmdphi);
imagTommy[0] = factor*(imagBarry[0]*cosmdphi-realBarry[0]*sinmdphi);
realTommy[1] = factor*(realBarry[1]*cosmdphi+imagBarry[1]*sinmdphi)*(r*r-2*r+a*a)/(r*r);
imagTommy[1] = factor*(imagBarry[1]*cosmdphi-realBarry[1]*sinmdphi)*(r*r-2*r+a*a)/(r*r);
realTommy[2] = factor*(realBarry[2]*cosmdphi+imagBarry[2]*sinmdphi)/r;
imagTommy[2] = factor*(imagBarry[2]*cosmdphi-realBarry[2]*sinmdphi)/r;
realTommy[3] = factor*(realBarry[3]*cosmdphi+imagBarry[3]*sinmdphi)/(r*sin(th));
imagTommy[3] = factor*(imagBarry[3]*cosmdphi-realBarry[3]*sinmdphi)/(r*sin(th));
realTommy[4] = factor*(realBarry[4]*cosmdphi+imagBarry[4]*sinmdphi)*(r*r-2*r+a*a)*(r*r-2*r+a*a)/(r*r*r*r);
imagTommy[4] = factor*(imagBarry[4]*cosmdphi-realBarry[4]*sinmdphi)*(r*r-2*r+a*a)*(r*r-2*r+a*a)/(r*r*r*r);
realTommy[5] = factor*(realBarry[5]*cosmdphi+imagBarry[5]*sinmdphi)*(r*r-2*r+a*a)/(r*r*r);
imagTommy[5] = factor*(imagBarry[5]*cosmdphi-realBarry[5]*sinmdphi)*(r*r-2*r+a*a)/(r*r*r);
realTommy[6] = factor*(realBarry[6]*cosmdphi+imagBarry[6]*sinmdphi)*(r*r-2*r+a*a)/(r*r*r*sin(th));
imagTommy[6] = factor*(imagBarry[6]*cosmdphi-realBarry[6]*sinmdphi)*(r*r-2*r+a*a)/(r*r*r*sin(th));
realTommy[7] = factor*(realBarry[7]*cosmdphi+imagBarry[7]*sinmdphi)/(r*r);
imagTommy[7] = factor*(imagBarry[7]*cosmdphi-realBarry[7]*sinmdphi)/(r*r);
realTommy[8] = factor*(realBarry[8]*cosmdphi+imagBarry[8]*sinmdphi)/(r*r*sin(th));
imagTommy[8] = factor*(imagBarry[8]*cosmdphi-realBarry[8]*sinmdphi)/(r*r*sin(th));
realTommy[9] = factor*(realBarry[9]*cosmdphi+imagBarry[9]*sinmdphi)/(r*r*sin(th)*sin(th));
imagTommy[9] = factor*(imagBarry[9]*cosmdphi-realBarry[9]*sinmdphi)/(r*r*sin(th)*sin(th));
}

int SeffTommy(int m, double a, double r, double th, double realBarry[10],  double imagBarry[10], double realTommy[10], double imagTommy[10]){
double factor = 1.0/(2*M_PI);
double rplus = 1+sqrt(1-a*a);
double rminus = 1-sqrt(1-a*a);
double mdphi = m*a/(rplus-rminus)*log((r-rplus)/(r-rminus));
double cosmdphi = cos(mdphi);
double sinmdphi = sin(mdphi);
realTommy[0] = factor*(realBarry[0]*cosmdphi+imagBarry[0]*sinmdphi)*(-((r*(a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/((a*a+r*r)*(a*a+r*r))));
imagTommy[0] = factor*(imagBarry[0]*cosmdphi-realBarry[0]*sinmdphi)*(-((r*(a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/((a*a+r*r)*(a*a+r*r))));
realTommy[1] = factor*(realBarry[1]*cosmdphi+imagBarry[1]*sinmdphi)*(-0.5*((a*a+(-2+r)*r)*(a*a+(-2+r)*r)*(a*a+2*r*r+a*a*cos(2*th)))/(r*(a*a+r*r)*(a*a+r*r)));
imagTommy[1] = factor*(imagBarry[1]*cosmdphi-realBarry[1]*sinmdphi)*(-0.5*((a*a+(-2+r)*r)*(a*a+(-2+r)*r)*(a*a+2*r*r+a*a*cos(2*th)))/(r*(a*a+r*r)*(a*a+r*r)));
realTommy[2] = factor*(realBarry[2]*cosmdphi+imagBarry[2]*sinmdphi)*(-(((a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/((a*a+r*r)*(a*a+r*r))));
imagTommy[2] = factor*(imagBarry[2]*cosmdphi-realBarry[2]*sinmdphi)*(-(((a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/((a*a+r*r)*(a*a+r*r))));
realTommy[3] = factor*(realBarry[3]*cosmdphi+imagBarry[3]*sinmdphi)*(-(((a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/(sin(th)*(a*a+r*r)*(a*a+r*r))));
imagTommy[3] = factor*(imagBarry[3]*cosmdphi-realBarry[3]*sinmdphi)*(-(((a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/(sin(th)*(a*a+r*r)*(a*a+r*r))));
realTommy[4] = factor*(realBarry[4]*cosmdphi+imagBarry[4]*sinmdphi)*(-0.5*((a*a+(-2+r)*r)*(a*a+(-2+r)*r)*(a*a+(-2+r)*r)*(a*a+2*r*r+a*a*cos(2*th)))/(r*r*r*(a*a+r*r)*(a*a+r*r)));
imagTommy[4] = factor*(imagBarry[4]*cosmdphi-realBarry[4]*sinmdphi)*(-0.5*((a*a+(-2+r)*r)*(a*a+(-2+r)*r)*(a*a+(-2+r)*r)*(a*a+2*r*r+a*a*cos(2*th)))/(r*r*r*(a*a+r*r)*(a*a+r*r)));
realTommy[5] = factor*(realBarry[5]*cosmdphi+imagBarry[5]*sinmdphi)*(-0.5*((a*a+(-2+r)*r)*(a*a+(-2+r)*r)*(a*a+2*r*r+a*a*cos(2*th)))/(r*r*(a*a+r*r)*(a*a+r*r)));
imagTommy[5] = factor*(imagBarry[5]*cosmdphi-realBarry[5]*sinmdphi)*(-0.5*((a*a+(-2+r)*r)*(a*a+(-2+r)*r)*(a*a+2*r*r+a*a*cos(2*th)))/(r*r*(a*a+r*r)*(a*a+r*r)));
realTommy[6] = factor*(realBarry[6]*cosmdphi+imagBarry[6]*sinmdphi)*(-(((a*a+(-2+r)*r)*(a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/(sin(th)*r*r*(a*a+r*r)*(a*a+r*r))));
imagTommy[6] = factor*(imagBarry[6]*cosmdphi-realBarry[6]*sinmdphi)*(-(((a*a+(-2+r)*r)*(a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/(sin(th)*r*r*(a*a+r*r)*(a*a+r*r))));
realTommy[7] = factor*(realBarry[7]*cosmdphi+imagBarry[7]*sinmdphi)*(-(((a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/(r*(a*a+r*r)*(a*a+r*r))));
imagTommy[7] = factor*(imagBarry[7]*cosmdphi-realBarry[7]*sinmdphi)*(-(((a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/(r*(a*a+r*r)*(a*a+r*r))));
realTommy[8] = factor*(realBarry[8]*cosmdphi+imagBarry[8]*sinmdphi)*(-(((a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/(sin(th)*r*(a*a+r*r)*(a*a+r*r))));
imagTommy[8] = factor*(imagBarry[8]*cosmdphi-realBarry[8]*sinmdphi)*(-(((a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/(sin(th)*r*(a*a+r*r)*(a*a+r*r))));
realTommy[9] = factor*(realBarry[9]*cosmdphi+imagBarry[9]*sinmdphi)*(-(((a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/(sin(th)*sin(th)*r*(a*a+r*r)*(a*a+r*r))));
imagTommy[9] = factor*(imagBarry[9]*cosmdphi-realBarry[9]*sinmdphi)*(-(((a*a+(-2+r)*r)*(r*r+a*a*cos(th)*cos(th)))/(sin(th)*sin(th)*r*(a*a+r*r)*(a*a+r*r))));
}