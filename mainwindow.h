#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QtCharts>
#include <QChartView>
#include <QtCharts/QLineSeries>
#include "QVector"
#include <math.h>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();


private slots:
    void on_DrawButton_clicked();

    void on_ClearButton_clicked();


private:
    Ui::MainWindow *ui;
    double X_Left = -1.0;
    double X_Right = 1.0;
    double mu_1 = 0.0;
    double mu_2 = 0.0;

    double delta = 1.0/9.0;

    size_t n_param = 9;
    size_t n_c_param = 9;

    double step = (X_Right - X_Left) / n_param;

    bool data_changed = true;
    bool cleared = true;
    char task = 0;
    char func = 0;

    double maxErr_fu = 0.0;
    double maxErr_fu_x = 0.0;
    double maxErr_fud1 = 0.0;
    double maxErr_fud1_x = 0.0;
    double maxErr_fud2 = 0.0;
    double maxErr_fud2_x = 0.0;

    QChartView *chartView = new QChartView();
    QLineSeries *series_main = new QLineSeries();
    QLineSeries *series_pnts = new QLineSeries();
    QLineSeries *series_err = new QLineSeries();
    QChart *chart = new QChart();
    QValueAxis *axisX = new QValueAxis();
    QValueAxis *axisY = new QValueAxis();

    // #novectorwednesday
    double* func_array = nullptr;
    double* step_array = nullptr;
    double* A_array = nullptr;
    double* B_array = nullptr;
    double* C_array = nullptr;
    double* FI_array = nullptr;
    double* c_coeff_alpha_array = nullptr;
    double* c_coeff_beta_array = nullptr;
    double* c_coeff_array = nullptr;
    double* a_coeff_array = nullptr;
    double* b_coeff_array = nullptr;
    double* d_coeff_array = nullptr;

    void Clear_Part();
    void Clear_Full();
    void GetData();

    void Fill_FuncNStep_Array();
    void Fill_C_Mat();
    void Solve();

    double TrueFunction(double point);
    double TrueFunctionD(double point);
    double TrueFunctionD2(double point);

    double Spline(double point);
    double SplineD(double point);
    double SplineD2(double point);
};
#endif // MAINWINDOW_H
