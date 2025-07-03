#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <iterator>
#include <algorithm>
#include <windows.h>

using namespace std;


class FileOpenException : public exception {
public:
    explicit FileOpenException(const string& filename) : msg("Не вдалося відкрити файл: " + filename) {}
    const char* what() const noexcept override { return msg.c_str(); }
private:
    string msg;
};

class InputException : public exception {
public:
    explicit InputException(const string& message) : msg("Помилка введення: " + message) {}
    const char* what() const noexcept override { return msg.c_str(); }
private:
    string msg;
};

class MathException : public exception {
public:
    explicit MathException(const string& message) : msg("Математична помилка: " + message) {}
    const char* what() const noexcept override { return msg.c_str(); }
private:
    string msg;
};

class DataRangeException : public exception {
public:
    explicit DataRangeException(const string& message) : msg("Дані поза діапазоном: " + message) {}
    const char* what() const noexcept override { return msg.c_str(); }
private:
    string msg;
};

bool isNaN(double x);
map<double, pair<double, double>> readData(const string& filename);
pair<double, double> interpolate(const map<double, pair<double, double>>& values, double x);
double T(double x, const map<double, pair<double, double>>& data);
double U(double x, const map<double, pair<double, double>>& data);
double fun1(double x, double y, double z, const map<double, pair<double, double>>& data);
double Grs(double x, double y, const map<double, pair<double, double>>& data);
double Srz(double x, double y, double z, const map<double, pair<double, double>>& data);
double Gold(double x, double y);
double Glr(double x, double y);
double fun2(double x, double y, double z, const map<double, pair<double, double>>& data);
double Grs1(double x, double y, const map<double, pair<double, double>>& data);
double Gold1(double x, double y);
double Glr1(double x, double y);
double fun3(double x, double y, double z);

bool isNaN(double x) {
    return x != x;
}

map<double, pair<double, double>> readData(const string& filename) {
    map<double, pair<double, double>> data;
    ifstream file(filename);
    if (!file.is_open()) {
        throw FileOpenException(filename);
    }
    double x, t, u;
    while (file >> x >> t >> u) {
        data[x] = { t, u };
    }
    file.close();
    if (data.empty()) {
        throw FileOpenException("Файл " + filename + " порожній або недійсний");
    }
    return data;
}

pair<double, double> interpolate(const map<double, pair<double, double>>& values, double x) {
    if (isNaN(x)) throw MathException("Недійсне введення (NaN)");
    if (abs(x) > 1e6) throw DataRangeException("Вихід за межі діапазону");

    auto it = values.lower_bound(x);
    if (it == values.end()) {
        --it;
        if (it == values.begin() && x < it->first) throw DataRangeException("Значення нижче діапазону");
        return it->second;
    }
    if (it == values.begin()) {
        if (x > it->first) throw DataRangeException("Значення вище діапазону");
        return it->second;
    }

    map<double, pair<double, double>>::const_iterator prevIt = it;
    if (prevIt != values.begin()) {
        --prevIt;
    }
    else {
        prevIt = it;
    }

    if (x == it->first) return it->second;
    if (x == prevIt->first) return prevIt->second;

    double x1 = prevIt->first, y1_t = prevIt->second.first, y1_u = prevIt->second.second;
    double x2 = it->first, y2_t = it->second.first, y2_u = it->second.second;
    double t = y1_t + (y2_t - y1_t) * (x - x1) / (x2 - x1);
    double u = y1_u + (y2_u - y1_u) * (x - x1) / (x2 - x1);
    return { t, u };
}

double T(double x, const map<double, pair<double, double>>& data) {
    return interpolate(data, x).first;
}

double U(double x, const map<double, pair<double, double>>& data) {
    return interpolate(data, x).second;
}

double fun1(double x, double y, double z, const map<double, pair<double, double>>& data) {
    if (isNaN(x) || isNaN(y) || isNaN(z)) throw MathException("Недійсне введення (NaN)");
    if (abs(x) > 1e6 || abs(y) > 1e6 || abs(z) > 1e6) throw DataRangeException("Вихід за межі діапазону");
    return x * Grs(y, z, data) + y * Grs(x, z, data) + 0.33 * y * Grs(x, z, data);
}

double Grs(double x, double y, const map<double, pair<double, double>>& data) {
    if (isNaN(x) || isNaN(y)) throw MathException("Недійсне введення (NaN)");
    if (abs(x) > 1e6 || abs(y) > 1e6) throw DataRangeException("Вихід за межі діапазону");
    return 0.1389 * Srz(x + y, Gold(x, y), Glr(x, x * y), data) +
        1.8389 * Srz(-y, Gold(y, x / 5), Glr(5 * x, x * y), data) +
        0.83 * Srz(x - 0.9, Glr(y, x / 5), Gold(5 * y, y), data);
}

double Srz(double x, double y, double z, const map<double, pair<double, double>>& data) {
    if (isNaN(x) || isNaN(y) || isNaN(z)) throw MathException("Недійсне введення (NaN)");
    if (abs(x) > 1e6 || abs(y) > 1e6 || abs(z) > 1e6) throw DataRangeException("Вихід за межі діапазону");
    return (T(x, data) + U(z, data) - T(y, data)) * x - y;
}

double Gold(double x, double y) {
    if (isNaN(x) || isNaN(y)) throw MathException("Недійсне введення (NaN)");
    if (x > y && y != 0) return 0.15;
    if (x <= y && x != 0) return 0.1;
    return 0;
}

double Glr(double x, double y) {
    if (isNaN(x) || isNaN(y)) throw MathException("Недійсне введення (NaN)");
    if (x < 1) return y;
    if (x >= 1 && y < 1) return y;
    if (sqrt(x * x + y * y) - 4 < 0.1) return 1;
    return sqrt(x * x + y * y) - 4;
}

double fun2(double x, double y, double z, const map<double, pair<double, double>>& data) {
    if (isNaN(x) || isNaN(y) || isNaN(z)) throw MathException("Недійсне введення (NaN)");
    if (abs(x) > 1e6 || abs(y) > 1e6 || abs(z) > 1e6) throw DataRangeException("Вихід за межі діапазону");
    return x * Grs1(x, y, data) + y * Grs1(y, z, data) + z * Grs1(z, x, data);
}

double Grs1(double x, double y, const map<double, pair<double, double>>& data) {
    if (isNaN(x) || isNaN(y)) throw MathException("Недійсне введення (NaN)");
    if (abs(x) > 1e6 || abs(y) > 1e6) throw DataRangeException("Вихід за межі діапазону");
    return 0.14 * Srz(x + y, Gold1(x, y), Glr1(x, x * y), data) +
        1.83 * Srz(-y, Gold1(y, x / 5), Glr1(4 * x, x * y), data) +
        0.83 * Srz(x, Glr1(y, x / 4), Gold1(4 * y, y), data);
}

double Gold1(double x, double y) {
    if (isNaN(x) || isNaN(y)) throw MathException("Недійсне введення (NaN)");
    if (x > y && y > 0.1) return 0.15;
    if (x <= y && x > 0.1) return 0.1;
    return 0;
}

double Glr1(double x, double y) {
    if (isNaN(x) || isNaN(y)) throw MathException("Недійсне введення (NaN)");
    if (x < 1) return y;
    if (y >= 1) return 1;
    return 0;
}

double fun3(double x, double y, double z) {
    if (isNaN(x) || isNaN(y) || isNaN(z)) throw MathException("Недійсне введення (NaN)");
    if (abs(x) > 1e6 || abs(y) > 1e6 || abs(z) > 1e6) throw DataRangeException("Вихід за межі діапазону");
    return 1.3498 * z + 2.2362 * y - 2.348 * x * y;
}

int main() {

    SetConsoleOutputCP(1251);
    SetConsoleCP(1251);

    double x, y, z;
    cout << "Введіть x, y, z: ";
    if (!(cin >> x >> y >> z)) {
        cout << "Помилка: Невірний формат введення" << endl;
        return 1;
    }

    map<double, pair<double, double>> data;

    try {
        string filename;
        if (abs(x) <= 1) {
            filename = "X_1_1.dat";
            if (x > 0) {
                filename = "X1_00.dat";
            }
        }
        else if (x < -1 || x > 1) {
            filename = "X00_1.dat";
        }
        else {
            throw DataRangeException("Немає відповідного файлу даних");
        }

        data = readData(filename);
        cout << "Використовуються дані з " << filename << endl;

        cout << "fun1(x, y, z) = " << fun1(x, y, z, data) << endl;
        cout << "fun2(x, y, z) = " << fun2(x, y, z, data) << endl;
        cout << "fun3(x, y, z) = " << fun3(x, y, z) << endl;
    }
    catch (const FileOpenException& e) {
        cout << "Помилка: " << e.what() << endl;
        return 1;
    }
    catch (const InputException& e) {
        cout << "Помилка: " << e.what() << endl;
        return 1;
    }
    catch (const MathException& e) {
        cout << "Помилка: " << e.what() << endl;
        return 1;
    }
    catch (const DataRangeException& e) {
        cout << "Помилка: " << e.what() << endl;
        return 1;
    }
    catch (const exception& e) {
        cout << "Помилка: Невідома помилка - " << e.what() << endl;
        return 1;
    }

    return 0;
}