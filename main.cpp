#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <ctime>

using namespace std;

void print_vector(vector<double> & v);

vector<double> TDMA(int N, const vector<double> & a, const vector<double> & c, const vector<double> & b, const vector<double> & f)
{

    vector<double> s(N+1, 0.0);
    vector<double> m(N+1, 0.0);
    vector<double> x(N+1);

    s[1] = -b[0]/c[0];
    m[1] = f[0]/c[0];
    cout << "* ";
    for (int i = 1; i <= N - 1; ++i) {
        cout << "* ";
        s[i+1] = -b[i]/(a[i]*s[i] + c[i]);
        m[i+1] = (f[i] - a[i]*m[i])/(a[i]*s[i] + c[i]);

    }
    x[N] = (f[N] - a[N]*m[N])/(a[N]*s[N] + c[N]);
    for (int i = N-1; i >= 0; i--) {
        x[i] = s[i+1]*x[i+1] + m[i+1];
    }
    return x;

}

vector<double> tridiagonal_solver(int N, vector<double> & a, vector<double> & b,
        vector<double> & c, vector<double> & f){

    cout << "TD beg \n";
    cout << setprecision(15);
    int n = N+1;
    vector<double> x(n);

    for(int i=1; i<n; i++){

        double m = a[i-1]/b[i-1];
        b[i] -= m*c[i-1];
        f[i] -= m*f[i-1];
    }
    cout << " b: \n";
    for(int i=1; i<n; i++){
        cout << b[i] << " ";
    }
    cout << "\n f: \n";
    for(int i=1; i<n; i++){
        cout << f[i] << " ";
    }
    // solve for last x value
    x[n-1] = f[n-1]/b[n-1];

    // solve for remaining x values by back substitution
    for(int i=n-2; i >= 0; i--)
        x[i] = (f[i]- c[i]*x[i+1])/b[i];

    cout << "TD end \n";
    return x;
}

void nonTR_solver(int N, vector<double> & a, vector<double> & b,
                                  vector<double> & c, vector<double> & f)
{
    vector<vector<double>> MAT(N, vector<double>(N));

    for (int i = 0; i < N+1; ++i) {
        for (int j = 0; j < N+1; ++j) {
            MAT[i][j] = 0;
        }
    }
    MAT[0][0] = b[0];
    MAT[0][1] = c[0];
    for (int i = 1; i < N; ++i) {
        MAT[i][i] = b[i];
        MAT[i][i-1] = a[i-1];
        MAT[i][i+1] = c[i];
    }
    MAT[N][N-1] = a[N-1];
    MAT[N][N] = b[N];
    for (int i = 0; i < N+1; ++i) {
        for (int j = 0; j < N+1; ++j) {
            cout << " " << MAT[j][i];
        }
        cout <<endl;
    }
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
    for (int k = 0; k <= 80; k++){ //while(abs(S_k)>eps_Sn){
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


void solve_imp(int N, int M, double F, double S)
{
    double h = 1.0/double(N);
    double tau = h*h*F;

    ofstream NUM("/Users/marinastruleva/Desktop/TEST/NUM_SOL.txt");
    ofstream EXACT("/Users/marinastruleva/Desktop/TEST/EXACT_SOL.txt");
    NUM << setprecision(8);
    EXACT << setprecision(8);

    double eps_Tastr = -1;

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
    int m = 0;
    //while (eps_Tastr > 0.00001){
    for (int m = 0; m < M; ++m) {
        eps_Tastr = -1;
        u[0] = 0;
        for (int j = 1; j < N; ++j) {
            u[j] = F*(u_prev[j+1] - 2*u_prev[j] + u_prev[j-1]) + u_prev[j];
        }
        u[N] = (2*h - u[N-2] + 4*u[N-1])/3;
        //cout << "u:";
        //print_vector(u);
        for (int k = 0; k < N+1; ++k) {
            eps_Tastr = max(eps_Tastr, abs(k*h-u[k]));
        }
        //cout << "eps(" <<  m << ") = " << eps_Tastr << endl;

        u_prev = u;

        for(int i=0; i<N+1; i++) {
            NUM << u[i] << "\t"; // << u[i] << "\n";
            EXACT << f_cond(i*h, m*tau) << "\t"; // << u[i] << "\n";
        }

        NUM << "\n";
        EXACT << "\n";
        m++;
    }
    //for (int m = 0; m < M; ++m) {

    cout << "\nu(x," << M*tau << "): ";
    print_vector(u);
}


double EPS_solve_imp(int N, int M, double F, double S)
{
    double h = 1.0/double(N);
    cout << "\n N = "<<N << endl;
    double tau = h*h*F;

    double eps_T = -1;

    vector<double> u(N+1), u_prev(N+1);

    for (int j = 0; j < N + 1; ++j) {
        u_prev[j] = f_bound(j*h);
    }
    cout << "\nu(x,0):";
    print_vector(u_prev);
    //while (eps_Tastr > 0.00001){
    for (int m = 0; m < M; ++m) {
        double eps = -1;
        u[0] = 0;
        for (int j = 1; j < N; ++j) {
            u[j] = F*(u_prev[j+1] - 2*u_prev[j] + u_prev[j-1]) + u_prev[j];
        }
        u[N] = (2*h - u[N-2] + 4*u[N-1])/3;
        u_prev = u;
        for (int i = 0; i < N+1; ++i) {
            eps = max(eps, abs(u[i] - f_cond(i*h, m*tau)));
        }
        eps_T = max(eps_T, eps);
        //cout << "eps_T = " << eps_T << endl;
    }
    //cout << "\nu(x," << M*tau << "): ";
    print_vector(u);
    return eps_T;
}

void solve_exp(int N, int M, double F, double S)
{
    /*
     * S = 1
     * F = tau/h^2
     * s
     *
     * */

    double h = 1.0/double(N);
    double tau = h*h*F;

    ofstream NUM("/Users/marinastruleva/Desktop/TEST/NUM_SOL.txt");
    ofstream EXACT("/Users/marinastruleva/Desktop/TEST/EXACT_SOL.txt");
    NUM << setprecision(8);
    EXACT << setprecision(8);


    double eps_Tastr = -1;

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
        DG[k] = 1 + 2*F;
        LW[k] = -F;
        UP[k] = -F;
    }
    DG[0] = 1;
    UP[0] = 0;
    LW[N] = -4 + (1+2*F)/F;
    DG[N+1] = 2;

    for (int m = 0; m < M; ++m) {
        RH[0] = 0;
        for (int j = 1; j < N; ++j) {
            RH[j] = u_prev[j];
        }
        RH[N] = 2*h + u_prev[N-1]/F;

        u = tridiagonal_solver(N, LW, DG, UP, RH);
        //cout << "u:";
        //print_vector(u);
        for (int k = 0; k < N+1; ++k) {
            eps_Tastr = max(eps_Tastr, abs(u_prev[k]-u[k]));
        }
        u_prev = u;

        cout << "eps(" <<  m << ") = " << eps_Tastr << endl;
        for(int i=0; i<N+1; i++) {
            NUM << u[i] << "\t"; // << u[i] << "\n";
            EXACT << f_cond(i*h, m*tau) << "\t"; // << u[i] << "\n";
        }

        NUM << "\n";
        EXACT << "\n";
    }
    cout << "\nu(x," << M*tau << "): ";
    print_vector(u);
}

int main()
{
    const clock_t begin_time = clock();
    cout.precision(10);
    cout << "\t\t ======= BEGIN =====\n";
    int N = 10;
    int M = 15;
    double F = .025;
    double S = 0;
    double h = 1/double(N);
    cout << "\nfor N = " << N << ", M =  " << M << ", F = " << F << ", S = " << S;

    double eps_T;
    double eps_T2;
    double eps_T4;

    eps_T = EPS_solve_imp(N, M, F, S);
    eps_T2 = EPS_solve_imp(N, M*2, F/2, S);
    eps_T4 = EPS_solve_imp(N, M*4, F/4, S);

    cout << "\neps_T = " << eps_T;
    cout << "\neps_T2 = " << eps_T2;
    cout << "\neps_T4 = " << eps_T4;

    double res_t = (eps_T - eps_T2)/(eps_T2 - eps_T4);
    cout << "\n\tT conv rate: " << res_t;

    double eps_H;
    double eps_H2;
    double eps_H4;

    eps_H = EPS_solve_imp(N, M, F, S);
    eps_H2 = EPS_solve_imp(N*2, M, F, S);
    eps_H4 = EPS_solve_imp(N*4, M, F, S);

    cout << "\neps_H = " << eps_H;
    cout << "\neps_H2 = " << eps_H2;
    cout << "\neps_H4 = " << eps_H4;

    double res_h1 = (eps_H - eps_H2);
    double res_h2 = (eps_H2 - eps_H4);

    cout << "\n\t eps_H - eps_H2: " << res_h1;
    cout << "\n\t eps_H2 - eps_H4: " << res_h2;
    cout << "\n\tH conv rate: " << res_h1/res_h2;



    cout << "\n\n TEST \n";

    vector<double> a = {5, 2, 1, 0};
    vector<double> b = {2, 10, 4, -1, 1};
    vector<double> c = {0, 4, 1, 1};
    vector<double> d = {1.3, 4.5, 7.8, 0, 2};
    vector<double> t;

    //nonTR_solver(4, a, b, c, d);

    cout << "\n" << b[0];
    t = tridiagonal_solver(a.size(), a, b, c, d);
    print_vector(t);
    cout << "\n";

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

