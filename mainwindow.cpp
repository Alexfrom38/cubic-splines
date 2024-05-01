#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <cmath>
#include <QDebug>
#include <QColor>
#include <fstream>
#include <iomanip>
#include <iostream>

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::Clear_Part() {
    series_main->detachAxis(axisX);
    series_main->detachAxis(axisY);
    series_pnts->detachAxis(axisX);
    series_pnts->detachAxis(axisY);
    series_err->detachAxis(axisX);
    series_err->detachAxis(axisY);
    chart->removeAxis(axisX);
    chart->removeAxis(axisY);
    ui->chartLayout->removeWidget(chartView);
    chart->removeSeries(series_main);
    chart->removeSeries(series_pnts);
    chart->removeSeries(series_err);
    disconnect(chartView);
    disconnect(chart);
    delete chart;
    delete chartView;
    delete axisX;
    delete axisY;

    chartView = new QChartView();
    chart = new QChart();
    axisX = new QValueAxis();
    axisY = new QValueAxis();
}

void MainWindow::Clear_Full() {
    series_main->detachAxis(axisX);
    series_main->detachAxis(axisY);
    series_pnts->detachAxis(axisX);
    series_pnts->detachAxis(axisY);
    series_err->detachAxis(axisX);
    series_err->detachAxis(axisY);
    chart->removeAxis(axisX);
    chart->removeAxis(axisY);
    ui->chartLayout->removeWidget(chartView);
    chart->removeSeries(series_main);
    chart->removeSeries(series_pnts);
    chart->removeSeries(series_err);
    disconnect(chartView);
    disconnect(chart);
    delete chart;
    delete chartView;
    delete axisX;
    delete axisY;

    if(series_main->count() != 0) series_main->clear();
    if(series_pnts->count() != 0) series_pnts->clear();
    if(series_err->count() != 0) series_err->clear();

    chartView = new QChartView();
    chart = new QChart();
    axisX = new QValueAxis();
    axisY = new QValueAxis();

    ui->aux_info->clear();

    ui->tableWidget->setRowCount(0);
    ui->tableWidget->setColumnCount(0);

    ui->tableWidget_sub->setRowCount(0);
    ui->tableWidget_sub->setColumnCount(0);

    if (A_array != nullptr) { delete[] A_array; }
    if (B_array != nullptr) { delete[] B_array; }
    if (C_array != nullptr) { delete[] C_array; }
    if (FI_array != nullptr) { delete[] FI_array; }
    if (c_coeff_alpha_array != nullptr) { delete[] c_coeff_alpha_array; }
    if (c_coeff_beta_array != nullptr) { delete[] c_coeff_beta_array; }
    if (c_coeff_array != nullptr) { delete[] c_coeff_array; }
    if (a_coeff_array != nullptr) { delete[] a_coeff_array; }
    if (b_coeff_array != nullptr) { delete[] b_coeff_array; }
    if (d_coeff_array != nullptr) { delete[] d_coeff_array; }

    A_array = nullptr;
    B_array = nullptr;
    C_array = nullptr;
    FI_array = nullptr;
    c_coeff_alpha_array = nullptr;
    c_coeff_beta_array = nullptr;
    c_coeff_array = nullptr;
    a_coeff_array = nullptr;
    b_coeff_array = nullptr;
    d_coeff_array = nullptr;

    data_changed = true;
    cleared = true;

    maxErr_fu = 0.0;
    maxErr_fu_x = 0.0;
    maxErr_fud1 = 0.0;
    maxErr_fud1_x = 0.0;
    maxErr_fud2 = 0.0;
    maxErr_fud2_x = 0.0;
}

void MainWindow::GetData() {
    size_t temp_i = 0;

    data_changed = false;

    temp_i = ui->IterationCount->text().toInt();
    if (temp_i != n_param) {
        n_param = (temp_i <= (size_t)0) ? (1) : (temp_i);
        data_changed = true;
        step = (X_Right - X_Left) / n_param;
    }

    temp_i = ui->IterationCount_c->text().toInt();
    if (temp_i != n_c_param) {
        n_c_param = (temp_i <= (size_t)1) ? (2) : (temp_i);
        data_changed = true;
        delta = (X_Right - X_Left) / n_c_param;
    }

    if(ui->mainTask->isChecked())
    {
        if(task != 0){
            task = 0;
            data_changed = true;
        }
    }
    else if(ui->testTask->isChecked())
    {
        if(task != 1){
            task = 1;
            data_changed = true;
        }
    }
    else if(ui->oscillTask->isChecked())
    {
        if(task != 2){
            task = 2;
            data_changed = true;
        }
    }
    else this->close();

    switch (task) {
    case 0:
        X_Left = 2.0;
        X_Right = 4.0;
        break;
    case 1:
        X_Left = -1.0;
        X_Right = 1.0;
        break;
    case 2:
        X_Left = 2.0;
        X_Right = 4.0;
        break;
    }

    if(ui->funct->isChecked())
    {
        if(func != 0){
            func = 0;
            data_changed = true;
        }
    }
    else if(ui->diriv->isChecked())
    {
        if(func != 1){
            func = 1;
            data_changed = true;
        }
    }
    else if(ui->diriv2->isChecked())
    {
        if(func != 2){
            func = 2;
            data_changed = true;
        }
    }
    else this->close();

    if(ui->naturalGU->isChecked())
    {
        mu_1 = 0.0;
        mu_2 = 0.0;
        data_changed = true;
    }
    else if(ui->accurateGU->isChecked())
    {
        mu_1 = TrueFunctionD2(X_Left);
        mu_2 = TrueFunctionD2(X_Right);
        data_changed = true;
    }
    else this->close();
}

void MainWindow::on_DrawButton_clicked() {
    QTableWidgetItem *tbl;

    GetData();
    if(data_changed) { Clear_Full(); }
    else { Clear_Part(); }

    ui->chartLayout->addWidget(chartView);
    chart->legend()->hide();
    if (data_changed || cleared) {
        Fill_FuncNStep_Array();
        Fill_C_Mat();
        Solve();

        ui->tableWidget->setColumnCount(7);
        ui->tableWidget->setRowCount(n_param+1);
        ui->tableWidget->setHorizontalHeaderLabels(QStringList()<<"i"<<"x[i-1]"<<"x[i]"<<"a"<<"b"<<"c"<<"d");
        for (size_t i = 0; i < n_param + 1; i++) {
            tbl = new QTableWidgetItem(QString::number(i));
            ui->tableWidget->setItem(i, 0, tbl);

            tbl = new QTableWidgetItem(QString::number(X_Left + ((i-1>0)?(i-1):(0)) * step_array[((i-1>0)?(i-1):(0))]));
            ui->tableWidget->setItem(i, 1, tbl);

            tbl = new QTableWidgetItem(QString::number(X_Left + i * step_array[i]));
            ui->tableWidget->setItem(i, 2, tbl);

            tbl = new QTableWidgetItem(QString::number(a_coeff_array[i]));
            ui->tableWidget->setItem(i, 3, tbl);

            tbl = new QTableWidgetItem(QString::number(b_coeff_array[i]));
            ui->tableWidget->setItem(i, 4, tbl);

            tbl = new QTableWidgetItem(QString::number(c_coeff_array[i]));
            ui->tableWidget->setItem(i, 5, tbl);

            tbl = new QTableWidgetItem(QString::number(d_coeff_array[i]));
            ui->tableWidget->setItem(i, 6, tbl);
        }

        ui->tableWidget_sub->setColumnCount(11);
        ui->tableWidget_sub->setRowCount((int)((X_Right - X_Left)/delta) + 1);
        ui->tableWidget_sub->setHorizontalHeaderLabels(QStringList()<<"j"<<"x[j]"<<"F(x[j])"<<"S(x[j])"<<"F(x[j])-S(x[j])"<<"F'(x[j])"<<"S'(x[j])"<<"F'(x[j]-S'(x[j]))"<<"F''(x[j])"<<"S''(x[j])"<<"F''(x[j]-S''(x[j]))");
        for (size_t j = 0; j < (size_t)((X_Right - X_Left)/delta) + 1; j++) {
            double pnt = X_Left + j * delta;
            double fu = Spline(pnt);
            double fud1 = SplineD(pnt);
            double fud2 = SplineD2(pnt);
            double tfu = TrueFunction(pnt);
            double tfud1 = TrueFunctionD(pnt);
            double tfud2 = TrueFunctionD2(pnt);

            tbl = new QTableWidgetItem(QString::number(j));
            ui->tableWidget_sub->setItem(j, 0, tbl);

            tbl = new QTableWidgetItem(QString::number(pnt));
            ui->tableWidget_sub->setItem(j, 1, tbl);

            tbl = new QTableWidgetItem(QString::number(tfu));
            ui->tableWidget_sub->setItem(j, 2, tbl);

            tbl = new QTableWidgetItem(QString::number(fu));
            ui->tableWidget_sub->setItem(j, 3, tbl);

            tbl = new QTableWidgetItem(QString::number(tfu - fu));
            ui->tableWidget_sub->setItem(j, 4, tbl);

            tbl = new QTableWidgetItem(QString::number(tfud1));
            ui->tableWidget_sub->setItem(j, 5, tbl);

            tbl = new QTableWidgetItem(QString::number(fud1));
            ui->tableWidget_sub->setItem(j, 6, tbl);

            tbl = new QTableWidgetItem(QString::number(tfud1 - fud1));
            ui->tableWidget_sub->setItem(j, 7, tbl);

            tbl = new QTableWidgetItem(QString::number(tfud2));
            ui->tableWidget_sub->setItem(j, 8, tbl);

            tbl = new QTableWidgetItem(QString::number(fud2));
            ui->tableWidget_sub->setItem(j, 9, tbl);

            tbl = new QTableWidgetItem(QString::number(tfud2 - fud2));
            ui->tableWidget_sub->setItem(j, 10, tbl);

            switch (func) {
            case 0:
                series_pnts->append(pnt, fu);
                series_err->append(pnt, std::abs(tfu - fu) * 5.0);
                break;
            case 1:
                series_pnts->append(pnt, fud1);
                series_err->append(pnt, std::abs(tfud1 - fud1) * 5.0);
                break;
            case 2:
                series_pnts->append(pnt, fud2);
                series_err->append(pnt, std::abs(tfud2 - fud2) * 5.0);
                break;
            }

            if (maxErr_fu < std::abs(tfu - fu)) {
                maxErr_fu = std::abs(tfu - fu);
                maxErr_fu_x = pnt;
            }
            if (maxErr_fud1 < std::abs(tfud1 - fud1)) {
                maxErr_fud1 = std::abs(tfud1 - fud1);
                maxErr_fud1_x = pnt;
            }
            if (maxErr_fud2 < std::abs(tfud2 - fud2)) {
                maxErr_fud2 = std::abs(tfud2 - fud2);
                maxErr_fud2_x = pnt;
            }
        }

        for (size_t i = 0; i < n_param + 1; i++) {
            switch (func) {
            case 0:
                series_main->append(X_Left + i * step_array[i], TrueFunction(X_Left + i * step_array[i]));
                break;
            case 1:
                series_main->append(X_Left + i * step_array[i], TrueFunctionD(X_Left + i * step_array[i]));
                break;
            case 2:
                series_main->append(X_Left + i * step_array[i], TrueFunctionD2(X_Left + i * step_array[i]));
                break;
            }
        }

        ui->tableWidget->verticalHeader()->setVisible(false);
        ui->tableWidget_sub->verticalHeader()->setVisible(false);

        ui->aux_info->append(QString("Узлов (осн. сетка): ") + QString::number(n_param + 1));
        ui->aux_info->append(QString("Узлов (контр. сетка): ") + QString::number(n_c_param + 1));
        ui->aux_info->append(QString("max|F(x[j])-S(x[j])| = ") + QString::number(maxErr_fu) + QString(" (при x = ") + QString::number(maxErr_fu_x) + QString(")"));
        ui->aux_info->append(QString("max|F'(x[j])-S'(x[j])| = ") + QString::number(maxErr_fud1) + QString(" (при x = ") + QString::number(maxErr_fud1_x) + QString(")"));
        ui->aux_info->append(QString("max|F''(x[j])-S''(x[j])| = ") + QString::number(maxErr_fud2) + QString(" (при x = ") + QString::number(maxErr_fud2_x) + QString(")"));
    }

    QPen pen_1(QRgb(0xff0000));
    pen_1.setWidth(1);
    series_main->setPen(pen_1);
    chart->addSeries(series_main);

    QPen pen_2(QRgb(0x00ff00));
    pen_2.setWidth(1);
    series_pnts->setPen(pen_2);
    chart->addSeries(series_pnts);

    QPen pen_3(QRgb(0x0000ff));
    pen_3.setWidth(1);
    series_err->setPen(pen_3);
    chart->addSeries(series_err);

    chart->setTitle("График");

    axisX->setTitleText("x");
    axisX->setLabelFormat("%f");
    axisX->setTickCount(11);
    axisX->setRange(X_Left, X_Right);
    chart->addAxis(axisX, Qt::AlignBottom);
    series_pnts->attachAxis(axisX);
    series_main->attachAxis(axisX);
    series_err->attachAxis(axisX);

    axisY->setTitleText("v");
    axisY->setLabelFormat("%f");
    axisY->setTickCount(11);
    //axisY->setRange(0.0, 10.0);
    chart->addAxis(axisY, Qt::AlignLeft);
    series_pnts->attachAxis(axisY);
    series_main->attachAxis(axisY);
    series_err->attachAxis(axisY);

    chartView->setChart(chart);

    cleared = false;
}

void MainWindow::on_ClearButton_clicked() {
    Clear_Full();
}

void MainWindow::Fill_FuncNStep_Array() {
    if (func_array != nullptr) { delete[] func_array; }
    if (step_array != nullptr) { delete[] step_array; }
    func_array = new double[n_param + 1];
    step_array = new double[n_param + 1];
    for(size_t i = 0; i < n_param + 1; i++) {
        double x = X_Left + step * i;
        func_array[i] = TrueFunction(x);
        step_array[i] = step;
    }
}

void MainWindow::Fill_C_Mat() {
    if (A_array != nullptr) { delete[] A_array; }
    if (B_array != nullptr) { delete[] B_array; }
    if (C_array != nullptr) { delete[] C_array; }
    if (FI_array != nullptr) { delete[] FI_array; }
    if (c_coeff_alpha_array != nullptr) { delete[] c_coeff_alpha_array; }
    if (c_coeff_beta_array != nullptr) { delete[] c_coeff_beta_array; }
    if (c_coeff_array != nullptr) { delete[] c_coeff_array; }
    if (a_coeff_array != nullptr) { delete[] a_coeff_array; }
    if (b_coeff_array != nullptr) { delete[] b_coeff_array; }
    if (d_coeff_array != nullptr) { delete[] d_coeff_array; }
    A_array = new double[n_param + 1];
    B_array = new double[n_param + 1];
    C_array = new double[n_param + 1];
    FI_array = new double[n_param + 1];
    c_coeff_alpha_array = new double[n_param + 1];
    c_coeff_beta_array = new double[n_param + 1];
    c_coeff_array = new double[n_param + 1];
    a_coeff_array = new double[n_param + 1];
    b_coeff_array = new double[n_param + 1];
    d_coeff_array = new double[n_param + 1];

    A_array[0] = 0.0;
    B_array[0] = 0.0;
    C_array[0] = 1.0;
    FI_array[0] = mu_1;

    A_array[n_param] = 0.0;
    B_array[n_param] = 0.0;
    C_array[n_param] = 1.0;
    FI_array[n_param] = mu_2;

    for (size_t i = 1; i < n_param; i++) {
        FI_array[i] = 6.0 * (((func_array[i + 1] - func_array[i]) / step_array[i + 1]) - ((func_array[i] - func_array[i - 1]) / step_array[i]));
    }

    for (size_t i = 1; i < n_param; i++) {
        A_array[i] = step_array[i];
        B_array[i] = step_array[i + 1];
        C_array[i] = 2.0 * (step_array[i] + step_array[i + 1]);
    }

    c_coeff_alpha_array[0] = 0.0;
    c_coeff_beta_array[0] = mu_1;
    c_coeff_array[0] = mu_1;

    c_coeff_alpha_array[n_param] = 0.0;
    c_coeff_beta_array[n_param] = mu_2;
    c_coeff_array[n_param] = mu_2;

    a_coeff_array[0] = 0.0;
    b_coeff_array[0] = 0.0;
    d_coeff_array[0] = 0.0;
}

void MainWindow::Solve() {
    for (size_t i = 1; i < n_param; i++){
        c_coeff_alpha_array[i] = -B_array[i] / (C_array[i] + A_array[i] * c_coeff_alpha_array[i - 1]);
        c_coeff_beta_array[i] = (FI_array[i] - A_array[i] * c_coeff_beta_array[i - 1]) / (C_array[i] + c_coeff_alpha_array[i - 1] * A_array[i]);
    }

    for (size_t i = n_param - 1; i > 0; i--){
        c_coeff_array[i] = c_coeff_alpha_array[i] * c_coeff_array[i + 1] + c_coeff_beta_array[i];
    }

    for (size_t i = 1; i < n_param + 1; i++){
        a_coeff_array[i] = func_array[i];
        d_coeff_array[i] = (c_coeff_array[i] - c_coeff_array[i - 1]) / step_array[i];
        b_coeff_array[i] = (func_array[i] - func_array[i - 1]) / step_array[i] + c_coeff_array[i] * step_array[i] / 3.0 + c_coeff_array[i - 1] * step_array[i] / 6.0;
    }
}

double MainWindow::Spline(double point) {
    size_t i = std::min((size_t)((point - X_Left) / step), n_param - 1);
    double _point = point - (X_Left + ((double)i + 1.0) * step);
    return (a_coeff_array[i + 1] + (b_coeff_array[i + 1] * _point) + (c_coeff_array[i + 1] / 2.0 * pow(_point, 2)) + (d_coeff_array[i + 1] / 6.0 * pow(_point, 3)));
}

double MainWindow::SplineD(double point) {
    size_t i = std::min((size_t)((point - X_Left) / step), n_param - 1);
    double _point = point - (X_Left + ((double)i + 1.0) * step);
    return (b_coeff_array[i + 1] + (c_coeff_array[i + 1] * _point) + (d_coeff_array[i + 1] / 2.0 * pow(_point, 2)));
}

double MainWindow::SplineD2(double point) {
    size_t i = std::min((size_t)((point - X_Left) / step), n_param - 1);
    double _point = point - (X_Left + ((double)i + 1.0) * step);
    return (c_coeff_array[i + 1] + (d_coeff_array[i + 1] * _point));
}

double MainWindow::TrueFunction(double point) {
    switch (task) {
    case 0:
        return (log(point + 1.0) / point);
    case 1:
        return ((point < 0.0) ? (point * point * point + 3.0 * point * point) : (-point * point * point + 3.0 * point * point));
    case 2:
        return (log(point + 1.0) / point + cos(10.0 * point));
    }
}

double MainWindow::TrueFunctionD(double point) {
    switch (task) {
    case 0:
        return (1.0 / (point * point + point) - (log(point + 1) / (point * point)));
    case 1:
        return ((point < 0.0) ? ((3.0 * point * point) + 6.0 * point) : ((-3.0 * point * point) + 6.0 * point));
    case 2:
        return (1.0 / (point * point + point) - (log(point + 1) / (point * point)) - 10.0 * sin(10.0 * point));
    }
}

double MainWindow::TrueFunctionD2(double point) {
    switch (task) {
    case 0:
        return (-(2.0 * point + 1.0)/(point * point * point * point + 2.0 * point * point * point + point * point) - 1.0 / (point * point * point + point * point) + 2.0 * log(point + 1) / (point * point * point));
    case 1:
        return ((point < 0.0) ? (6.0 * point + 6.0) : (-6.0 * point + 6.0));
    case 2:
        return (-(2.0 * point + 1.0)/(point * point * point * point + 2.0 * point * point * point + point * point) - 1.0 / (point * point * point + point * point) + 2.0 * log(point + 1) / (point * point * point) - 100.0 * cos(10.0 * point));
    }
}
