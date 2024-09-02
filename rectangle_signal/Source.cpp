#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <complex>
using namespace std;

#define N_SECONDS 100

double rectangle_impulse(double x, double T, double t, double x_m, double phi) {
    x = x - phi;
    if (x < 0)
        while (x < 0) x = x + T;
    while (x > T) x = x - T;
    if (x <= t) return x_m;
    else return 0;
}

void ComplexBitReverse(complex<double>* data, int size) {
    int middle = size / 2, j = 0;
    for (int i = 0; i < size - 1; ++i) {
        if (i < j) swap(data[i], data[j]); // меняем элементы местами
        int k = middle;
        while (k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

void FftDit(complex<double>* data, int size, int sizeLog2, int dir) {
    ComplexBitReverse(data, size);
    int ptsInLeftDft, ptsInRightDft = 1;
    for (int stage = 1; stage <= sizeLog2; ++stage) {
        ptsInLeftDft = ptsInRightDft;   // установить ptsInLeftDFT = 2**(stage-1)
        ptsInRightDft *= 2;             // установить ptsInRightDFT = 2**stage
        complex<double> twiddle = complex<double>(1.0, 0.0); // поворачивающий множ.
        double trigArg = M_PI / ptsInLeftDft;
        // dir == 1 для прямого преобразования, dir == -1 для обратного
        complex<double> wFactor = complex<double>(cos(trigArg), -sin(trigArg) * dir);
        for (int butterflyPos = 0; butterflyPos < ptsInLeftDft; ++butterflyPos) {
            for (int topNode = butterflyPos; topNode < size; topNode += ptsInRightDft) {
                int botNode = topNode + ptsInLeftDft;
                complex<double> temp = data[botNode] * twiddle;
                data[botNode] = data[topNode] - temp;
                data[topNode] += temp;
            }  // конец цикла по topNode

            twiddle *= wFactor;
        } // конец цикла "бабочка"
    } // конец цикла stage
}

void fill_csv_rectangle(double F, double quantization_max, double T, double t, double phi) {
    int N_SAMPLES = N_SECONDS * F;
    double quantization_min = 0.;
    uint8_t quantization_levels_num = 64;
    uint8_t* digital_signal = new uint8_t[N_SAMPLES];
    for (int i = 0; i < N_SAMPLES; i++) digital_signal[i] = (round(rectangle_impulse(i / F, T, t, quantization_max, phi) * (quantization_levels_num - 1))) /
        (quantization_max - quantization_min);
    const char csv_file_name[64] = "data_0.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    for (size_t i = 0; i < N_SAMPLES; ++i)
    {
        double signal_val = double(digital_signal[i]) * (quantization_max - quantization_min) /
            (quantization_levels_num - 1) + quantization_min;

        csv_file << signal_val * pow(-1, i) << "\n";
    }
    csv_file.close();
    delete[] digital_signal;
}


int main() {
    double F = 40.; // частота дискретизации
    double quantization_max = 1.;
    double T = 6.;
    double t = 4.;
    double phi = 0.;
    int size = pow(2, round(log(N_SECONDS * F + 1) / log(2)));

    fill_csv_rectangle(F, quantization_max, T, t, phi);

    complex<double>* data = new complex<double>[size];
    std::ifstream input_data;
    input_data.open("data_0.csv");

    string line;
    int i = 0;
    while (getline(input_data, line)) {
        double re = stod(line);
        complex<double> buf(re, 0.);
        data[i] = buf;
        i++;
    }
    if (i < size) {
        for (int j = i; j < size; j++) data[j] = complex<double>(0, 0);
    }

    input_data.close();

    FftDit(data, size, log(size) / log(2), 1);
    for (int i = 0; i < size; i++) {
        data[i] = data[i] / (complex<double>)size;
        //cout << data[i] << endl;
    }

    std::ofstream res;
    res.open("data.csv");
    res << "time, amplitude of signal" << "\n";
    for (int i = size / 2; i < size; i++) {
        res << (i - size / 2) * F / size << "," << sqrt(data[i].real() * data[i].real() + data[i].imag() * data[i].imag()) << "\n";
    }

    res.close();
}