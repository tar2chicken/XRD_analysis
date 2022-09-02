#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
using namespace std;




// Data of each peak
// データ型の定義
struct peakdata {
    int number;
    double theta2;
    double angle;
    double sin4;
    double sd;
    int si;
    double a;
    double nr;
};





int main() {
    // Specify iofile name
    // ファイル名の設定
    string ifname, rfname, pfname;
    std::cout << "Enter the name of the input file." << endl;
    cin >> ifname;
    std::cout << endl;
    rfname = ifname + "_result";
    pfname = ifname + "_plot";
    ifstream ifile(ifname);
    ofstream rfile(rfname);
    ofstream pfile(pfname);




    // Declaration
    // 宣言
    vector<peakdata> peaks;
    const double pi = acos(-1);




    // Confirm the wavelength
    // 波長を確認
    std::cout << "The default wavelength of incident X-rays is 1.5405 \u00C5." << endl;
    string yn;
    double lambda;
    for (int i = 0; i < 10; i++) {
        std::cout << "Do you want to continue? [Y/n]";
        cin >> yn;
        if (yn == "Y" || yn == "y" || yn == "YES" || yn == "Yes" || yn == "yes") {
            lambda = 1.5405;
            break;
        } else if (yn == "n" || yn == "N" || yn == "no" || yn == "No" || yn == "NO") {
            std::cout << "Enter wavelength in \u00C5." << endl;
            cin >> lambda;
            if (lambda > 0.0) {
                break;
            } else {
                std::cout << "Invalid value." << endl;
            }
        }
        if (i < 9) {
            std::cout << "Sorry, try again." << endl << endl;
        } else {
            std::cout << "Set the value to 1.5405 \u00C5." << endl;
            lambda = 1.5405;
        }
    }
    std::cout << endl;





    // Read data from ifile
    // データの読み込み
    int n;
    ifile >> n;
    for (int i = 0; i < n; i++) {
        peakdata ipeak = {0, 0, 0, 0, 0, 0, 0, 0};
        ifile >> ipeak.number >> ipeak.theta2;
        ipeak.angle = ipeak.theta2 * pi / 360.0;
        ipeak.sin4 = 4.0 * pow(sin(ipeak.angle), 2) / pow(lambda, 2);
        ipeak.nr = (1/sin(ipeak.angle) + 1/ipeak.angle) * pow(cos(ipeak.angle), 2);
        peaks.push_back(ipeak);
    }





    // Determine the lattice type and parameter
    // 格子の種類の決定
    string lattice_t;
    double lattice_p;
    double ratio21;
    ratio21 = peaks.at(1).sin4 / peaks.at(0).sin4;
    if (ratio21 > 2.33) {
        lattice_t = "diamond structure";
        lattice_p = 3.0;
    } else if (ratio21 < 1.67) {
        lattice_t = "FCC (face-centered cubic)";
        lattice_p = 3.0;
    } else if (n >= 7) {
        double ratio71;
        ratio71 = peaks.at(6).sin4 / peaks.at(0).sin4;
        if ( ratio71 < 7.5) {
                lattice_t = "BCC (body-centered cubic)";
                lattice_p = 2.0;
        } else {
            lattice_t = "SC (simple cubic)";
            lattice_p = 1.0;
        }
    } else {
        lattice_t = "uncertain";
        lattice_p = 1.0;
    }





    // Determine the index and the lattice constant of each peak
    // ミラー指数と格子定数の決定
    for (int i = 0; i < n; i++) {
        peaks.at(i).sd = peaks.at(i).sin4 * lattice_p / peaks.at(0).sin4;
        peaks.at(i).si = round(peaks.at(i).sd);
        peaks.at(i).a = sqrt(peaks.at(i).si / peaks.at(i).sin4);
    }






    // Refine the lattice constant
    // 格子定数の精密決定
    double F, F2, FA, A, A2;
    F = 0.0;
    F2 = 0.0;
    FA = 0.0;
    A = 0.0;
    A2 = 0.0;
    for (int i = 0; i < n; i++) {
        F = F + peaks.at(i).nr;
        F2 = F2 + pow(peaks.at(i).nr, 2);
        FA = FA + (peaks.at(i).nr * peaks.at(i).a);
        A = A + peaks.at(i).a;
        A2 = A2 + pow(peaks.at(i).a, 2);
    }
    double slope, interc, sigma0, sigmas, sigmai, sigma0s, sigma0i;
    slope = (n*FA - F*A) / (n*F2 - pow(F, 2));
    interc = (F2*A - FA*F) / (n*F2 - pow(F, 2));
    sigma0 = sqrt((A2 - n*pow(interc, 2) - F2*pow(slope, 2) - 2*slope*interc*F) / (n-2));
    sigmas = sqrt(n / (n*F2 - pow(F, 2)));
    sigmai = sqrt(F2 / (n*F2 - pow(F, 2)));
    sigma0s = sigma0 * sigmas;
    sigma0i = sigma0 * sigmai;





    // Write data to rfile
    // データの書き込み
    rfile << "No.  \u03B8(rad)  1/d^2(\u00C5^-2)  index  a(\u00C5)" << endl;
    rfile << "-------------------------------------" << endl;
    for (int i = 0; i < n; i++) {
        rfile << peaks.at(i).number << "  " << peaks.at(i).angle << "  " << peaks.at(i).sin4 << "  " << peaks.at(i).si << "  " << peaks.at(i).a << endl;
    }
    rfile << endl << "regression line  :  y = " << slope << "*x + " << interc << endl;
    rfile << "slope            :  " << slope << " \u00B1 " << sigma0s << endl;
    rfile << "intercept        :  " << interc << " \u00B1 " << sigma0i << endl;
    rfile << endl << endl << "The lattice type is " << lattice_t << "." << endl;
    rfile << "The lattice constant is " << interc << "\u00B1" << sigma0i << " \u00C5." << endl;
    rfile.close();






    // Write data to pfile
    // データの書き込み
    for (int i = 0; i < n; i++){
        pfile << peaks.at(i).nr << " " << peaks.at(i).a << endl;
    }
    pfile.close();





    // Finish
    // 終了
    std::cout << "Done!" << endl;
    std::cout << "The result is output to " << rfname << "." << endl;
    std::cout << "Data for refinement of the lattice constant is output to " << pfname << "." << endl;

    return 0;
}