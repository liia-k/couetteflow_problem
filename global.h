#ifndef GLOBAL_H
#define GLOBAL_H
#include <QVector>
#include <QtMath>

constexpr double Nav(6.02214129e23); // Avogadro constant, 1/mol
constexpr double UniversalGasConstant = 8.3144598; // R, J/(K*mol)
constexpr double kB=1.38064852e-23;    // Boltzmann constant, J/K
constexpr double kBE = 8.617333e-5;       // const Bolz elVolt/K

constexpr double hPlank = 6.62559e-34; // Planck constant, kg*m^2/s (Length^2Mass/Time)
constexpr double clight = 2.99792458e8; // m/s
constexpr double hc = hPlank*clight; //  kg*m^3/s^2  (L^3M/T^2=J*L)

struct solverParams
{
    int NumCell     = 0;     // Число расчетных ячеек
    double t_fin    = 0;     // Время выхода из решения
    double Ma       = 0;     // Число маха
    double Gamma    = 0;     // Показатель адиабаты
    double CFL      = 0;     // Число Куранта U*deltaT/deltaH<1
    double lambda   = 0;     // Длина свободного пробега
    int lambdaSol   = 0;     // Кол-во длин пробега для расчета
    int PlotIter    = 0;     // Кол-во итераций, через которое отрисовывается график
    int MaxIter     = 10000000; // максимальное кол-во шагов по времени
    int typePlot    = 0;     // Тип отображения на графике :
                             // 0 - абс. давление
                             // 1 - абс. плотность (общая)
                             // 2 - абс. скорость
                             // 3 - абс. температура
                             // 4 - абс. плотность компонента O2
                             // 5 - абс. плотность компонента O
    int typeRightBorder = 1; //  Тип граничного условия справа
                             //  1 -- первая модель
                             //  2 -- новая модель
    int typeLeftBorder = 1;  //  Тип граничного условия слева


    double h_length = 0.05;  // Расстояние между пластинами, м
    double temp_wall1 = 300; // Температура нижней пластины, К
    double temp_wall2 = 300; // Температура верхней пластины, К
    double vel_wall2 = 300;  // Скорость верхней пластины, м/с
    double accom_coeff1 = 0.2;// Коэффициент аккоммодации на нижней стенке
    double accom_coeff2 = 0.2;// Коэффициент аккоммодации на верхней стенке
    double rec_coeff1 = 0.02; // Коэффициент рекомбинации на нижней стенке
    double rec_coeff2 = 0.02; // Коэффициент рекомбинации на верхней стенке

};
class Matrix
{
private:
    QVector<double> data;
public:
    Matrix (QVector<double> val)
    {
        data = val;
    }
    Matrix()
    {

    }
    void clear()
    {
        data.clear();
    }
    Matrix (int len, double val = 0)
    {
        data = QVector<double> (len, val);
    }
    QVector<double>::iterator begin()
    {
        return data.begin();
    }
    QVector<double>::iterator end()
    {
        return data.end();
    }
    int size()
    {
        return data.size();
    }
    operator QVector<double>() const
    {
        return data;
    }
    double& operator [](int i)
    {
        static double elseValue;
        if(i < data.size())
            return data[i];
        else
            return elseValue;
    }
    double first()
    {
        return data.first();
    }
    double last()
    {
        return data.last();
    }
    void push_back(double val)
    {
        data.push_back(val);
    }
    void push_front(double val)
    {
        data.push_front(val);
    }
    void removeLast()
    {
        data.removeLast();
    }
    void removeFirst()
    {
        data.removeFirst();
    }
    void resize(int i)
    {
        data.resize(i);
    }
    Matrix operator /(QVector<double> div)
    {
        QVector<double> output;
        if(data.size() != div.size())
            return data;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] / div[i]);
        return output;
    }
    Matrix operator *(QVector<double> div)
    {
        QVector<double> output;
        if(data.size() != div.size())
            return data;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] * div[i]);
        return output;
    }
    Matrix operator +(QVector<double> div)
    {
        QVector<double> output;
        if(data.size() != div.size())
            return data;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] + div[i]);
        return output;
    }
    Matrix operator -(QVector<double> div)
    {
        QVector<double> output;
        if(data.size() != div.size())
            return data;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] - div[i]);
        return output;
    }
    Matrix operator /(double div)
    {
        QVector<double> output;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] / div);
        return output;
    }
    Matrix operator *(double div)
    {
        QVector<double> output;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] * div);
        return output;
    }
    Matrix operator +(double div)
    {
        QVector<double> output;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] + div);
        return output;
    }
    Matrix operator -(double div)
    {
        QVector<double> output;
        for(int i = 0; i < data.size(); i ++)
            output.push_back(data[i] - div);
        return output;
    }
    static Matrix POW(QVector<double> div, double param)
    {
        QVector<double> output;
        for(int i = 0; i < div.size(); i ++)
             output.push_back(pow(div[i],param));
        return output;
    }
    static Matrix SQRT(QVector<double> div)
    {
        QVector<double> output;
        for(int i = 0; i < div.size(); i ++)
             output.push_back(sqrt(div[i]));
        return output;
    }
    static Matrix REVERSE(QVector<double> div)
    {
        QVector<double> output;
        for(int i = 0; i < div.size(); i ++)
             output.push_back(1.0/div[i]);
        return output;
    }
};

#endif // GLOBAL_H
