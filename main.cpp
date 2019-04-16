#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <ctime>

using namespace std;

double eps_Sn = pow(10, -8);
double eps_X = -1;

void print_vector(vector<double> & v);

vector<double> tridiagonal_solver(int N, vector<double> & a, vector<double> & b, vector<double> & c, vector<double> & f);

double f_bound(double x){
    return sin(M_PI*x) + x;
}

vector<double> slava_TDMA(int n, vector<double> & A, vector<double> & C, vector<double> & B, vector<double> & F)
{

    int N = n;
    vector<double> s(N), m(N), y(N);

    s[1] = - B[0] / C[0];
    m[1] = F[0] / C[0];

    for (int i = 1; i <= N -1; i++)
    {
        s[i+1] = - B[i] / (A[i]*s[i] + C[i]);
        m[i+1] = (F[i] - A[i] * m[i])/(A[i] * s[i] + C[i]);
    }

    y[N] = (F[N] - A[N] * m[N])/(A[N] * s[N] + C[N]);

    for (int i = N -1; i >= 0; i--)
    {
        y[i] = s[i+1] * y[i+1] + m[i+1];
    }
    return y;
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
// основная функция
void solve(int N, int M, double tau, double SIGMA, bool show_here){

    double h = 1/double(N);
    double F = tau/h/h; //отношение tau/h^2 = F
    double S = SIGMA; // вес схемы


    // запись в файл. NUM - численное решение, EXACT - точное. по строкам
    ofstream NUM("/Users/marinastruleva/Desktop/TEST/NUM_SOL.txt");
    ofstream EXACT("/Users/marinastruleva/Desktop/TEST/EXACT_SOL.txt");
    NUM << setprecision(8);
    EXACT << setprecision(8);

    cout << " for N = " << N << " M = "<< M << ": \n\n";

    //F = tau/h/h;

    cout << " h = " << h << ", tau = " << tau << ", M*tau = " << M*tau << endl;

    vector<double> u(N+1), u_prev(N+1), u_real(N+1), u_lim(N+1);

    for (int i = 0; i < N+1; ++i) {
        u_prev[i] = f_bound(i*h);
        u_lim[i] = i*h;
        u_real[i] = f_cond(i*h, 0);
    }
    for(int i=0; i < N+1; i++) {
        NUM << i*h << "\t";// << u[i] << "\n";
        EXACT << i*h << "\t";// << u[i] << "\n";
        //file << endl; //<< u[i] << endl;
    }
    NUM << "\n";
    EXACT << "\n";

    // если включено show_here, то решение выводится на экран
    if (!show_here) {
        for (int i = 0; i < N+1; ++i) cout << i*h << "\t";
        cout << endl;
        cout << "u_prev \n";
        print_vector(u_prev);

        cout << "u_real \n";
        print_vector(u_real);

    }
    //print_vector(u_real);

//Initializing the TDMA
    //TDM_init(F, S, N);
    /* ------------ MAIN PART ------------ */

    vector<double> RH(N+1);

    vector<double> LW(N), DG(N+1), UP(N);

    DG[0] = 1;
    for(int i=1; i<N; i++){
        DG[i] = 1 + 2*F*S;
    }
    DG[N] = 2; //1
    LW[0] = -F*S;
    UP[0] = 0;
    for(int i=1; i<N-1; i++){
        LW[i] = -F*S;
        UP[i] = -F*S;
    }
    LW[N-1] = -4 + (1+2*F*S)/(F*S);
    UP[N-1] = -F*S;
    cout << " F = " << F << endl;
    if(!show_here) {
        cout << "\n ******* \n";
        cout << "LW: " << LW.size() << endl;
        print_vector(LW);
        cout << "DG: " << DG.size() << endl;
        print_vector(DG);
        cout << "UP: "<< UP.size()  << endl;
        print_vector(UP);
        cout << "\n ******* \n";
    }


    cout << "BEGIN " << F*S << endl;
//RHS initialization
// заполнение правой части RH и решение системы.
    for (int m = 0; m < M; ++m) { //потом будет условный цикл
        RH[0] = 0;
        for (int j = 1; j < N; j++) {
            RH[j] = u_prev[j] + (1-S)*F*(u_prev[j-1] - 2*u_prev[j] + u_prev[j+1]);
        }
        RH[N] = 2*h + RH[N-1]/S/F; // h;
        if (!show_here) {
            cout << endl;
            for (int i = 0; i < N+1; ++i) {
                cout << "RH[" << i << "]:" << RH[i] << endl;
            }

            //print_vector(RH);
        }
        u = tridiagonal_solver(N, LW, DG, UP, RH);
        if (!show_here) {
            cout << "u_prev: \n";
            print_vector(u_prev);

            cout << "u: \n";
            print_vector(u);
        }

        for(int i=0; i<N+1; i++) {
            NUM << u[i] << "\t"; // << u[i] << "\n";
            EXACT << f_cond(i*h, m*tau) << "\t"; // << u[i] << "\n";

        }

        NUM << "\n";
        EXACT << "\n";

        // Evaluating the errors

        u_prev = u;
        if (show_here){
            cout << "**\n*";
            cout << "u_prev: \n";
            print_vector(u_prev);

            cout << "u: \n";
            print_vector(u);
        }

        for(int i=0; i<N+1; i++) {
            eps_X = max(eps_X, abs(u_prev[i] - u_real[i]));
        }

        //u.clear();
        //RH.clear();
        cout << "m = " << m << endl;
    }

    cout << "!!" << endl;
// Writing results in file
    //file << "HI!";

    NUM.close();
    EXACT.close();


    cout << "*************" <<endl;
    // проверка согласованности размеров
    cout << "x: " << N+1 << endl;
    cout << "u: " << u_prev.size() << endl;
    cout << "LW: " << LW.size() << endl;
    cout << "DG: " << DG.size() << endl;
    cout << "UP: " << UP.size() << endl;

    RH.clear();
    u.clear();
/*
    string line;
    ifstream fin("/Users/marinastruleva/Desktop/TEST/NUM_SOL.txt");
    if(fin.is_open())
    {
        while ( getline (fin,line) )
        {
            cout << line << '\n';
        }
        fin.close();
    }
    else cerr<<"Unable to open file";
*/

}

int main(){
    const clock_t begin_time = clock();

    cout << setprecision(15);
    int N = 10;
    int M = 2;
    double h = 1.0/double(N);
    double tau = .25*h*h;
    double Sigma = 1;
    solve(N, M, tau, Sigma, false);

    // проверка прогонки. работает
    vector<double> a = {5, 2, 1, 0};
    vector<double> b = {2, 10, 4, -1, 1};
    vector<double> c = {0, 4, 1, 1};
    vector<double> d = {1.3, 4.5, 7.8, 0, 2};
    vector<double> t;

    cout << "\n" << b[0];
    t = slava_TDMA(a.size(), a, b, c, d);
    print_vector(t);

    cout << "\n time used = "<< float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    cout << "\n!END!";
    return 0;
}

// Printing vectors
void print_vector(vector<double> & v){
    cout<<endl;
    for(int i=0; i<v.size(); i++){
        cout << setprecision(7) << setw(5) << v[i] << "\t";
    }
    cout<<endl;
}

// Solving tridiagonal matrices


vector<double> tridiagonal_solver(int N, vector<double> & a, vector<double> & b, vector<double> & c, vector<double> & f){

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

