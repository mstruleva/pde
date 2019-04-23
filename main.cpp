#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <ctime>

using namespace std;

void print_vector(vector<double> & v);

vector<double> TDMA(int N, const vector<double> & lw, const vector<double> & dg, const vector<double> & up, const vector<double> & rh)
{

    vector<double> ac = lw, bc = dg, cc = up, dc = rh;

    int nf = N+1;
    vector<double> xc(nf);

    for (int i = 1; i < nf; ++i) {
        //cout << "\n " << ac[i-1] << " " << bc[i-1] << endl;
        double mc = ac[i-1]/bc[i-1];
        //cout << " mc: " << mc;
        bc[i] -= mc*cc[i-1];
        dc[i] -= mc*dc[i-1];
    }

    xc = bc;
    xc[nf-1] = dc[nf-1]/bc[nf-1];
    for (int i = nf-2; i >=0 ; i--) {
        xc[i] = (dc[i] - cc[i]*xc[i+1])/bc[i];
    }
    return xc;
}

double f_bound(double x)
{
    return sin(M_PI*x) + x;
}

double f_cond(double x, double t){
    double u = x;
    int k = 0;
    double l_k, A_k, S_k;
    l_k = pow((.5*M_PI + M_PI*k), 2); //
    A_k = (-8*cos(M_PI*k))/M_PI/(4*k*k + 4*k - 3); //
    S_k = A_k*sin((M_PI*k + .5*M_PI)*x)*exp(-l_k*t);
    u += S_k;
    //cout << S_k << endl;
    // k - количество слагаемых из ряда
    for (int k = 0; k <= 200; k++){ //while(abs(S_k)>eps_Sn){
        k +=1;
        l_k = pow((.5*M_PI + M_PI*k), 2);
        A_k = (-8*cos(M_PI*k))/M_PI/(4*k*k + 4*k - 3);
        S_k = A_k*sin((M_PI*k + .5*M_PI)*x)*exp(-l_k*t);
        u += S_k;
        //cout << k << ' ';
    }
    //cout << k << ' ';
    return  u;
}

double f_condX(double x, double t){
    return exp(-2.25*M_PI*M_PI*t) * sin(1.5*M_PI*x) + x;
}

void solve_imp(int N, int M, double F, double S)
{
    cout << "\nIMPLICIT:\n";
    double h = 1.0/double(N);
    double tau = h*h*F;

    ofstream NUM("/Users/marinastruleva/Desktop/TEST/NUM_SOL.txt");
    ofstream EXACT("/Users/marinastruleva/Desktop/TEST/EXACT_SOL.txt");
    ofstream EXACT_X("/Users/marinastruleva/Desktop/TEST/EXACT_X.txt");

    NUM << setprecision(8);
    EXACT << setprecision(8);

    //double eps_Tastr = -1;

    for(int i=0; i < N+1; i++) {
        NUM << i*h << "\t";// << u[i] << "\n";
        EXACT << i*h << "\t";// << u[i] << "\n";
        EXACT_X << i*h << "\t";// << u[i] << "\n";

    }
    NUM << "\n";
    EXACT << "\n";
    EXACT_X << "\n";


    vector<double> u(N+1), u_prev(N+1);

    for (int j = 0; j < N + 1; ++j) {
        u_prev[j] = f_bound(j*h);
    }
    cout << "\nu(x,0):";
    print_vector(u_prev);

    int m = 0;
    double eps = 1e-5;
    double eps_Tastr = 2;

    while (eps_Tastr > eps){
        eps_Tastr = -1;
        u[0] = 0;
        for (int j = 1; j < N; ++j) {
            u[j] = F*(u_prev[j+1] - 2*u_prev[j] + u_prev[j-1]) + u_prev[j];
        }
        u[N] = (2*h - u[N-2] + 4*u[N-1])/3;
        //cout << "u:";
        //print_vector(u);
        for (int k = 0; k < N+1; ++k) {
            eps_Tastr = max(eps_Tastr, abs(u_prev[k]-u[k]));
        }
        //cout << "eps(" <<  m << ") = " << eps_Tastr << endl;

        u_prev = u;

        for(int i=0; i<N+1; i++) {
            NUM << u[i] << "\t"; // << u[i] << "\n";
            EXACT << f_cond(i*h, m*tau) << "\t"; // << u[i] << "\n";
            EXACT_X << f_condX(i*h, m*tau) << "\t"; // << u[i] << "\n";
        }

        NUM << "\n";
        EXACT << "\n";
        EXACT_X << "\n";

        m++;
    }
    cout << "m = " << m;
    //for (int m = 0; m < M; ++m) {

    cout << "\nu(x," << M*tau << "): ";
    print_vector(u);
}


double EPS_solve_imp(int N, int M, double F, double S)
{
    double h = 1.0/double(N);
    //cout << "\n N = "<<N << endl;
    double tau = h*h*F;

    vector<double> u(N+1), u_prev(N+1);

    for (int j = 0; j < N + 1; ++j) {
        u_prev[j] = f_bound(j*h);
    }
    //cout << "\nu(x,0):";
    //print_vector(u_prev);
    double eps_stable = 1e-5;
    double eps_Tastr = 2;

    double eps_T = -1;

    int m = 0;
    while (eps_Tastr > eps_stable){
        eps_Tastr = -1;
        double eps = -1;
        u[0] = 0;
        for (int j = 1; j < N; ++j) {
            u[j] = F*(u_prev[j+1] - 2*u_prev[j] + u_prev[j-1]) + u_prev[j];
        }
        u[N] = (2*h - u[N-2] + 4*u[N-1])/3; //h + u[N-1];
        for (int i = 0; i < N+1; ++i) {
            eps = max(eps, abs(u[i] - f_cond(i*h, m*tau)));
        }
        for (int k = 0; k < N+1; ++k) {
            eps_Tastr = max(eps_Tastr, abs(u_prev[k]-u[k]));
        }
        //cout << "T* = " << eps_Tastr << endl;
        eps_T = max(eps_T, eps);
        m++;
        u_prev = u;
        //cout << "eps_T = " << eps_T << endl;
    }

    cout << "\nm = " << m <<endl;
    return eps_T;
}

double EPS_solve_exp(int N, int M, double F, double S)
{
    double h = 1.0/double(N);
    //cout << "\n N = "<<N << endl;
    double tau = h*h*F;

    double eps_T = -1;

    vector<double> u(N+1), u_prev(N+1);

    for (int j = 0; j < N + 1; ++j) {
        u_prev[j] = f_bound(j*h);
    }

    vector<double> LW(N), UP(N), DG(N+1), RH(N+1);
    for (int k = 0; k < N; ++k) {
        DG[k] = 1 + 2*F*S;
        LW[k] = -F*S;
        UP[k] = -F*S;
    }
    DG[0] = 1;
    UP[0] = 0;
    LW[N-1] = -4 + (1 + 2*F*S)/F/S; //0;
    DG[N] = 2; //1;

    for (int m = 0; m < M; ++m) {
        double eps = -1;
        RH[0] = 0;
        for (int j = 1; j < N; ++j) {
            RH[j] = u_prev[j] + (1 - S) * F * (u_prev[j - 1] - 2 * u_prev[j] + u_prev[j + 1]);
        }
        RH[N] = 2 * h + RH[N - 1] / F / S; //1;

        u = TDMA(N, LW, DG, UP, RH);
        u_prev = u;
        for (int i = 0; i < N + 1; ++i) {
            eps = max(eps, abs(u[i] - f_cond(i * h, m * tau)));
        }
        eps_T = max(eps_T, eps);
    }
    return eps_T;
}

void solve_exp(int N, int M, double F, double S)
{
    cout << "\nEXPLICIT:\n";
    double h = 1.0/double(N);
    double tau = h*h*F;

    ofstream NUM("/Users/marinastruleva/Desktop/TEST/NUM_SOL.txt");
    ofstream EXACT("/Users/marinastruleva/Desktop/TEST/EXACT_SOL.txt");
    NUM << setprecision(8);
    EXACT << setprecision(8);


    for(int i=0; i < N+1; i++) {
        NUM << i*h << "\t";// << u[i] << "\n";
        EXACT << i*h << "\t";// << u[i] << "\n";
    }
    NUM << "\n";
    EXACT << "\n";

    vector<double> u(N+1), u_prev(N+1);

    for (int j = 0; j < N + 1; ++j) {
        u_prev[j] = f_bound(j*h);
    }
    cout << "\nu(x,0):";
    print_vector(u_prev);

    vector<double> LW(N), UP(N), DG(N+1), RH(N+1);
    for (int k = 0; k < N; ++k) {
        DG[k] = 1 + 2*F*S;
        LW[k] = -F*S;
        UP[k] = -F*S;
    }
    DG[0] = 1;
    UP[0] = 0;
    LW[N-1] = -4 + (1 + 2*F*S)/F/S; //0;
    DG[N] = 2; //1;

    int m = 0;
    double eps = 1e-5;
    double eps_Tastr = 2;

    //for (int m = 0; m < M; ++m) {
    while (eps_Tastr > eps){
        eps_Tastr = -1;
        RH[0] = 0;
        for (int j = 1; j < N; ++j) {
            RH[j] = u_prev[j] + (1-S)*F*(u_prev[j-1] - 2*u_prev[j] + u_prev[j+1]);
        }
        RH[N] = 2*h + RH[N-1]/F/S; //1;

        u = TDMA(N, LW, DG, UP, RH);

        //cout << "u:";
        //print_vector(u);
        for (int k = 0; k < N+1; ++k) {
            eps_Tastr = max(eps_Tastr, abs(u_prev[k]-u[k]));
        }
        u_prev = u;

        //cout << "eps(" <<  m << ") = " << eps_Tastr << endl;
        for(int i=0; i<N+1; i++) {
            NUM << u[i] << "\t"; // << u[i] << "\n";
            EXACT << f_cond(i*h, m*tau) << "\t"; // << u[i] << "\n";
        }

        NUM << "\n";
        EXACT << "\n";
        //RH.clear();
        //u.clear();
        m++;
    }
    cout << "\nu(x," << m*tau << "): ";
    cout << "\n M = " << m << "\n tau*M = " << tau*m;

    print_vector(u);

    //
}

int main() {
    const clock_t begin_time = clock();
    cout << "\t\t ======= BEGIN =====\n";
    int N = 40;
    int M = 20;
    double F = .4;
    double S = .5;

    cout << "\nfor N = " << N << " M = " << M << " F = " << F << " S = " << S << endl;
    solve_imp(N, M, F, S);

    //eps count
    int IT = 12;
/*
    cout << "\nTau convergence rate from M = " << M << " to " <<M*pow(2, IT+1) << ":" << endl;
    for (int it = 2; it < IT; ++it) {
        double eps_T = EPS_solve_exp(N, int(M*pow(2,it-1)), F/(pow(2,it-1)), S);
        double eps_T2 = EPS_solve_exp(N, int(M*pow(2,it)), F/(pow(2,it)), S);
        double eps_T4 = EPS_solve_exp(N, int(M*pow(2,it+1)), F/(pow(2,it+1)), S);
        //cout << "\nM: " << M*pow(2,it-1) << "\t" << M*pow(2,it) << "\t"<< M*pow(2,it+1) << "\n";
        //cout << "F: " << F/(pow(2,it-1)) << "\t" << F/(pow(2,it)) << "\t"<< F/(pow(2,it+1)) << "\n";
        double res_t = (eps_T - eps_T2)/(eps_T2 - eps_T4);
        cout << "\n\t" << res_t << endl;
    }
    */
    /*
    cout << "\nH convergence rate from N = " << N << " to " <<N*pow(2, IT) << ":" << endl;
    double F_t = 1e-5;

    double eps_prev = EPS_solve_imp(N, M, F, S);
    for (int it = 1; it < IT + 1; ++it) {
        double eps_cur = EPS_solve_imp(int(N*pow(2, it)), M, F, S);
        double res = eps_cur/eps_prev;
        cout << N*pow(2, it) << "  \t" << res << endl;
        eps_prev = eps_cur;
    }
//    */
    cout << "\n\t\t ======= END =======";
    cout << "\n time used = "<< float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    return 0;
}

void print_vector(vector<double> & v) {
    cout << endl;
    for (int i = 0; i < v.size(); i++) {
        cout << setprecision(7) << setw(5) << v[i] << "\t";
    }
    cout << endl;
}