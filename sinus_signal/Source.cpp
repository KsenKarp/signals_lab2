#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <complex>
using namespace std;

#define N_SECONDS 100

void ComplexBitReverse(complex<double>* data, int size) {
    int middle = size / 2, j = 0;
    for (int i = 0; i < size - 1; ++i) {
        if (i < j) swap(data[i], data[j]); // ������ �������� �������
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
        ptsInLeftDft = ptsInRightDft;   // ���������� ptsInLeftDFT = 2**(stage-1)
        ptsInRightDft *= 2;             // ���������� ptsInRightDFT = 2**stage
        complex<double> twiddle = complex<double>(1.0, 0.0); // �������������� ����.
        double trigArg = M_PI / ptsInLeftDft;
        // dir == 1 ��� ������� ��������������, dir == -1 ��� ���������
        complex<double> wFactor = complex<double>(cos(trigArg), -sin(trigArg) * dir);
        for (int butterflyPos = 0; butterflyPos < ptsInLeftDft; ++butterflyPos) {
            for (int topNode = butterflyPos; topNode < size; topNode += ptsInRightDft) {
                int botNode = topNode + ptsInLeftDft;
                complex<double> temp = data[botNode] * twiddle;
                data[botNode] = data[topNode] - temp;
                data[topNode] += temp;
            }  // ����� ����� �� topNode

            twiddle *= wFactor;
        } // ����� ����� "�������"
    } // ����� ����� stage
}

void fill_csv_sin(double F, double quantization_min, double quantization_max, double phi, double T) {
    uint8_t quantization_levels_num = 64; // ���������� ������� ������������ �����������
    int N_SAMPLES = N_SECONDS * F + 1;
    uint8_t* digital_signal = new uint8_t[N_SAMPLES]; // ���� ����������� ��������� �������
    double k = M_PI * 2 / T; //(M_PI*2/T) - ����������� k � ������
    for (int i = 0; i < N_SAMPLES; i++) digital_signal[i] = (round((sin(i * k / F + phi) + 1) * (quantization_levels_num - 1))) /
        (quantization_max - quantization_min);
    const char csv_file_name[64] = "data_0.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    for (size_t i = 0; i < N_SAMPLES; ++i)
    {
        double signal_val = double(digital_signal[i]) * (quantization_max - quantization_min) /
            (quantization_levels_num - 1) + quantization_min;

        csv_file << signal_val * pow(-1, i) << "\n"; //����� � �������������� ����� ��� �� ���������
    }
    csv_file.close();
    delete[] digital_signal;
}


int main() {
    double F = 40.; // ������� �������������
    double quantization_min = -1., quantization_max = 1.;
    double phi = M_PI;
    double T = 2 * M_PI;
    int size = pow(2, round(log(N_SECONDS * F + 1) / log(2)));
    
    fill_csv_sin(F, quantization_min, quantization_max, phi, T);

    complex<double>* data = new complex<double>[size];
    std::ifstream input_data;
    input_data.open("data_0.csv");

    string line;
    int i = 0;
    while (getline(input_data, line)) {
        double re = stod(line);
        complex<double> buf (re, 0.);
        data[i] = buf;
        i++;
    }
    if (i < size) {
        for (int j = i; j < size; j++) data[j] = complex<double>(0, 0);
    }

    input_data.close();
    
    FftDit(data, size, log(size)/log(2), 1);
    for (int i = 0; i < size; i++) {
        data[i] = data[i] / (complex<double>)size;
        //cout << data[i] << endl;
    }

    std::ofstream res;
    res.open("data.csv");
    res << "time, amplitude of signal" << "\n";
    for (int i = size/2; i < size; i++) {
        res << (i - size / 2 ) * F / size << "," << sqrt(data[i].real() * data[i].real() + data[i].imag() * data[i].imag()) << "\n";
    }

    res.close();
}