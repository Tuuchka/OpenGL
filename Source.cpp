// Лабораторная работа №2
// М8О-305Б-20 Тучков Николай
// Вар. 22;  6 - гранная прямая усеченная пирамида
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <GL/glut.h>


void display();
void specialKeys();

const double pi = 3.141592;
double rotate_y = 0;
double rotate_x = 0;
double rotate_z = 0;
double transf = 1;
double R = 0.5; // радиус вписанной в 8-угольник окружность
double k = 1 + sqrt(2);
double t = 2 * R / k; // сторона 8-угольника
//double r = 0.3;
double sm = 0.1;
bool OZ = 0;
bool OZ2 = 0;
bool OX = 0;
bool OY = 0;

double VECT[3] = { 0, 0, 0 };

double VEC1[3];
double VEC2[3];

double FaceC[3];
double SCAL;
double T[4][4];

double K[3] = { 0, 0, 2 };

double X0[3];
double X1[3];
double X2[3];
double X3[3];
double X4[3];
double X5[3];
double X6[3];
double X7[3];
double X8[3];
double X9[3];
double X10[3];
double X11[3];
double X12[3];
double X13[3];
double X14[3];
double X15[3];

// Функция очистки переменных
void Clean() {
    for (int i = 0; i < 3; i++) {
        VEC1[i] = 0;
        VEC2[i] = 0;
        VECT[i] = 0;
        FaceC[i] = 0;
    }
    SCAL = 0;
}

double CUBE[16][4] = { -R - sm, t / 2, -R, 1,
                       -R - t * cos(3 * pi / 4) - sm, t * sin(3 * pi / 4) + t / 2, -R, 1,
                       -R - sm, -t / 2, -R, 1,
                       -R - t * cos(3 * pi / 4) - sm, -t * sin(3 * pi / 4) - t / 2, -R, 1,
                        R - sm, t / 2, -R, 1,
                        R + t * cos(3 * pi / 4) - sm, t * sin(3 * pi / 4) + t / 2, -R, 1,
                        R - sm, -t / 2, -R, 1,
                        R + t * cos(3 * pi / 4) - sm, -t * sin(3 * pi / 4) - t / 2, -R, 1,

                       -R + sm, t / 2, R, 1,
                       -R - t * cos(3 * pi / 4) + sm, t * sin(3 * pi / 4) + t / 2, R, 1,
                       -R + sm, -t / 2, R, 1,
                       -R - t * cos(3 * pi / 4) + sm, -t * sin(3 * pi / 4) - t / 2, R, 1,
                        R + sm, t / 2, R, 1,
                        R + t * cos(3 * pi / 4) + sm, t * sin(3 * pi / 4) + t / 2, R, 1,
                        R + sm, -t / 2, R, 1,
                        R + t * cos(3 * pi / 4) + sm, -t * sin(3 * pi / 4) - t / 2, R, 1 };

// Функция возвращения фигуры в изначальное положение
void CubeBack() {
    double A2[16][4] = { -R - sm, t / 2, -R, 1,
                       -R - t * cos(3 * pi / 4) - sm, t * sin(3 * pi / 4) + t / 2, -R, 1,
                       -R - sm, -t / 2, -R, 1,
                       -R - t * cos(3 * pi / 4) - sm, -t * sin(3 * pi / 4) - t / 2, -R, 1,
                        R - sm, t / 2, -R, 1,
                        R + t * cos(3 * pi / 4) - sm, t * sin(3 * pi / 4) + t / 2, -R, 1,
                        R - sm, -t / 2, -R, 1,
                        R + t * cos(3 * pi / 4) - sm, -t * sin(3 * pi / 4) - t / 2, -R, 1,

                       -R + sm, t / 2, R, 1,
                       -R - t * cos(3 * pi / 4) + sm, t * sin(3 * pi / 4) + t / 2, R, 1,
                       -R + sm, -t / 2, R, 1,
                       -R - t * cos(3 * pi / 4) + sm, -t * sin(3 * pi / 4) - t / 2, R, 1,
                        R + sm, t / 2, R, 1,
                        R + t * cos(3 * pi / 4) + sm, t * sin(3 * pi / 4) + t / 2, R, 1,
                        R + sm, -t / 2, R, 1,
                        R + t * cos(3 * pi / 4) + sm, -t * sin(3 * pi / 4) - t / 2, R, 1 };
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            CUBE[i][j] = A2[i][j];
            A2[i][j] = 0;
        }
    }
}

// Очистка матрицы Т
void ClearT() {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            T[i][j] = 0;
        }
    }
}

// Преобразование матрицы точек фигуры в векторы
void CUBEonVEC() {
    for (int i = 0; i < 3; i++) {
        X0[i] = CUBE[0][i];
        X1[i] = CUBE[1][i];
        X2[i] = CUBE[2][i];
        X3[i] = CUBE[3][i];
        X4[i] = CUBE[4][i];
        X5[i] = CUBE[5][i];
        X6[i] = CUBE[6][i];
        X7[i] = CUBE[7][i];
        X8[i] = CUBE[8][i];
        X9[i] = CUBE[9][i];
        X10[i] = CUBE[10][i];
        X11[i] = CUBE[11][i];
        X12[i] = CUBE[12][i];
        X13[i] = CUBE[13][i];
        X14[i] = CUBE[14][i];
        X15[i] = CUBE[15][i];
    }
}

// Функция векторного умножения
void VECxVEC(double VECT1[], double VECT2[]) {
    VECT[0] = VECT1[1] * VECT2[2] - VECT2[2] * VECT1[1];
    VECT[1] = VECT1[2] * VECT2[0] - VECT1[0] * VECT2[2];
    VECT[2] = VECT1[0] * VECT2[1] - VECT1[1] * VECT2[0];
}

// Функция скалярного умножения
void VECsVEC(double VECT1[], double VECT2[]) {
    SCAL = VECT1[0] * VECT2[0] + VECT1[1] * VECT2[1] + VECT1[2] * VECT2[2];
}

//Функция определения видимости граней для всех, кроме 7-ой и 8-ой грани
bool FACEisReal(double Y1[], double Y2[], double Y3[], double Y4[]) {
    VEC1[0] = Y1[0] - Y2[0];
    VEC1[1] = Y1[1] - Y2[1];
    VEC1[2] = Y1[2] - Y2[2];

    VEC2[0] = Y2[0] - Y3[0];
    VEC2[1] = Y2[1] - Y3[1];
    VEC2[2] = Y2[2] - Y3[2];

    VECxVEC(VEC1, VEC2); //VECT=

    for (int i = 0; i < 3; i++) {
        FaceC[i] = (Y1[i] + Y2[i] + Y3[i] + Y4[i]) / 4;
    }
    VECsVEC(FaceC, VECT); //SCAL=

    if (SCAL <= 0) {
        VECT[0] = -1 * VECT[0];
        VECT[1] = -1 * VECT[1];
        VECT[2] = -1 * VECT[2];
    }

    VECsVEC(K, VECT);
    if (SCAL < 0)
        return 0;
    else
        return 1;
}

//Функция определения видимости граней для 7-ой и 8-ой грани
bool FACEisReal2(double Y1[], double Y2[], double Y3[], double Y4[], double Y5[], double Y6[], double Y7[], double Y8[]) {
    VEC1[0] = Y1[0] - Y2[0];
    VEC1[1] = Y1[1] - Y2[1];
    VEC1[2] = Y1[2] - Y2[2];

    VEC2[0] = Y2[0] - Y3[0];
    VEC2[1] = Y2[1] - Y3[1];
    VEC2[2] = Y2[2] - Y3[2];

    VECxVEC(VEC1, VEC2);

    for (int i = 0; i < 3; i++) {
        FaceC[i] = (Y1[i] + Y2[i] + Y3[i] + Y4[i] + Y5[i] + Y6[i] + Y7[i] + Y8[i]) / 8;
    }
    VECsVEC(FaceC, VECT);

    if (SCAL < 0) {
        VECT[0] = -1 * VECT[0];
        VECT[1] = -1 * VECT[1];
        VECT[2] = -1 * VECT[2];
    }

    VECsVEC(K, VECT);
    if (SCAL <= 0)
        return 0;
    else
        return 1;
}
// Поворот фигуры вокруг оси X
void RotateX(int phi) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            T[i][j] = 0;
        }
    }
    double phiR = phi * 0.0175;
    T[0][0] = 1;
    T[1][1] = cos(phiR);
    T[1][2] = sin(phiR);
    T[2][1] = -sin(phiR);
    T[2][2] = cos(phiR);
    T[3][3] = 1;

    double A2[16][4];
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            double sum = 0;
            for (int k = 0; k < 4; k++) {
                sum += CUBE[i][k] * T[k][j];
            }
            A2[i][j] = sum;
        }
    }

    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            CUBE[i][j] = A2[i][j];
            A2[i][j] = 0;
        }
    }
    rotate_x = 0;
}
// Поворот фигуры вокруг оси Y
void RotateY(int teta) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            T[i][j] = 0;
        }
    }
    double tetaR = teta * 0.0175;
    T[0][0] = cos(tetaR);
    T[2][2] = cos(tetaR);
    T[2][0] = sin(tetaR);
    T[0][2] = -sin(tetaR);
    T[1][1] = 1;
    T[3][3] = 1;

    double A2[16][4];
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            double sum = 0;
            for (int k = 0; k < 4; k++) {
                sum += CUBE[i][k] * T[k][j];
            }
            A2[i][j] = sum;
        }
    }

    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            CUBE[i][j] = A2[i][j];
            A2[i][j] = 0;
        }
    }
    rotate_y = 0;
}

// Поворот фигуры вокруг оси Z
void RotateZ(int psi) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            T[i][j] = 0;
        }
    }
    double psiR = psi * 0.0175;
    T[0][0] = cos(psiR);
    T[1][1] = cos(psiR);
    T[0][1] = sin(psiR);
    T[1][0] = -sin(psiR);
    T[2][2] = 1;
    T[3][3] = 1;

    double A2[16][4];
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            double sum = 0;
            for (int k = 0; k < 4; k++) {
                sum += CUBE[i][k] * T[k][j];
            }
            A2[i][j] = sum;
        }
    }

    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            CUBE[i][j] = A2[i][j];
            A2[i][j] = 0;
        }
    }
    rotate_z = 0;
}
// Изменение размеров фигуры
void Transform(double flag) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            T[i][j] = 0;
        }
    }
    T[0][0] = flag;
    T[1][1] = flag;
    T[2][2] = flag;
    T[3][3] = flag;

    double A2[16][4];
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            double sum = 0;
            for (int k = 0; k < 4; k++) {
                sum += CUBE[i][k] * T[k][j];
            }
            A2[i][j] = sum;
        }
    }

    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            CUBE[i][j] = A2[i][j];
            A2[i][j] = 0;
        }
    }
    transf = 1;
}
//Отрисовка фигуры в изначальном положении
void OrtoZ() {
    CubeBack();
    OZ = 0;
}
//Отрисовка ортогональной проекции фигуры по оси Z
void OrtoZ2() {
    CubeBack();
    T[0][0] = 1;
    T[1][1] = cos(pi);
    T[1][2] = sin(pi);
    T[2][1] = -sin(pi);
    T[2][2] = cos(pi);
    T[3][3] = 1;

    double A2[16][4];
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            double sum = 0;
            for (int k = 0; k < 4; k++) {
                sum += CUBE[i][k] * T[k][j];
            }
            A2[i][j] = sum;
        }
    }

    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            CUBE[i][j] = A2[i][j];
            A2[i][j] = 0;
        }
    }
    OZ2 = 0;
}
// Отрисовка ортогональной проекции по оси X
void OrtoX() {
    CubeBack();
    T[0][0] = 1;
    T[1][1] = cos(-pi / 2);
    T[1][2] = sin(-pi / 2);
    T[2][1] = -sin(-pi / 2);
    T[2][2] = cos(-pi / 2);
    T[3][3] = 1;

    double A2[16][4];
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            double sum = 0;
            for (int k = 0; k < 4; k++) {
                sum += CUBE[i][k] * T[k][j];
            }
            A2[i][j] = sum;
        }
    }

    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 4; j++) {
            CUBE[i][j] = A2[i][j];
            A2[i][j] = 0;
        }
    }
    OX = 0;
}


//Отрисовка фигуры
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //glLoadIdentity();
    RotateX(rotate_x);
    RotateY(rotate_y);
    RotateZ(rotate_z);
    Transform(transf);

    if (OZ == 1)
        OrtoZ();
    else if (OZ2 == 1)
        OrtoZ2();
    else if (OX == 1)
        OrtoX();
    CUBEonVEC();
    //1 - ая грань
    if (FACEisReal(X8, X2, X0, X10) == 1) {
        glBegin(GL_LINES);
        glVertex3f(CUBE[0][0], CUBE[0][1], CUBE[0][2]);
        glVertex3f(CUBE[2][0], CUBE[2][1], CUBE[2][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[0][0], CUBE[0][1], CUBE[0][2]);
        glVertex3f(CUBE[8][0], CUBE[8][1], CUBE[8][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[2][0], CUBE[2][1], CUBE[2][2]);
        glVertex3f(CUBE[10][0], CUBE[10][1], CUBE[10][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[8][0], CUBE[8][1], CUBE[10][2]);
        glVertex3f(CUBE[10][0], CUBE[10][1], CUBE[10][2]);
        glEnd();
    }
    Clean();


    // 2-ая грань
    if (FACEisReal(X0, X1, X9, X8) == 1) {
        glBegin(GL_LINES);
        glVertex3f(CUBE[0][0], CUBE[0][1], CUBE[0][2]);
        glVertex3f(CUBE[1][0], CUBE[1][1], CUBE[1][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[1][0], CUBE[1][1], CUBE[1][2]);
        glVertex3f(CUBE[9][0], CUBE[9][1], CUBE[9][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[0][0], CUBE[0][1], CUBE[0][2]);
        glVertex3f(CUBE[8][0], CUBE[8][1], CUBE[8][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[8][0], CUBE[8][1], CUBE[8][2]);
        glVertex3f(CUBE[9][0], CUBE[9][1], CUBE[9][2]);
        glEnd();
    }
    Clean();

    // 3-ая грань
    if (FACEisReal(X1, X5, X13, X9) == 1) {
        glBegin(GL_LINES);
        glVertex3f(CUBE[1][0], CUBE[1][1], CUBE[1][2]);
        glVertex3f(CUBE[5][0], CUBE[5][1], CUBE[5][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[5][0], CUBE[5][1], CUBE[5][2]);
        glVertex3f(CUBE[13][0], CUBE[13][1], CUBE[13][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[1][0], CUBE[1][1], CUBE[1][2]);
        glVertex3f(CUBE[9][0], CUBE[9][1], CUBE[9][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[13][0], CUBE[13][1], CUBE[13][2]);
        glVertex3f(CUBE[9][0], CUBE[9][1], CUBE[9][2]);
        glEnd();
    }
    Clean();

    //4-ая грань
    if (FACEisReal(X4, X5, X13, X12) == 1) {
        glBegin(GL_LINES);
        glVertex3f(CUBE[5][0], CUBE[5][1], CUBE[5][2]);
        glVertex3f(CUBE[4][0], CUBE[4][1], CUBE[4][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[5][0], CUBE[5][1], CUBE[5][2]);
        glVertex3f(CUBE[13][0], CUBE[13][1], CUBE[13][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[4][0], CUBE[4][1], CUBE[4][2]);
        glVertex3f(CUBE[12][0], CUBE[12][1], CUBE[12][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[12][0], CUBE[12][1], CUBE[12][2]);
        glVertex3f(CUBE[13][0], CUBE[13][1], CUBE[13][2]);
        glEnd();
    }
    Clean();

    //5-ая грань
    if (FACEisReal(X4, X6, X14, X12) == 1) {
        glBegin(GL_LINES);
        glVertex3f(CUBE[4][0], CUBE[4][1], CUBE[4][2]);
        glVertex3f(CUBE[6][0], CUBE[6][1], CUBE[6][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[4][0], CUBE[4][1], CUBE[4][2]);
        glVertex3f(CUBE[12][0], CUBE[12][1], CUBE[12][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[6][0], CUBE[6][1], CUBE[6][2]);
        glVertex3f(CUBE[14][0], CUBE[14][1], CUBE[14][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[12][0], CUBE[12][1], CUBE[12][2]);
        glVertex3f(CUBE[14][0], CUBE[14][1], CUBE[14][2]);
        glEnd();
    }
        Clean();

    //6 - ая грань
    if (FACEisReal(X6, X7, X15, X14) == 1) {
        glBegin(GL_LINES);
        glVertex3f(CUBE[6][0], CUBE[6][1], CUBE[6][2]);
        glVertex3f(CUBE[7][0], CUBE[7][1], CUBE[7][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[6][0], CUBE[6][1], CUBE[6][2]);
        glVertex3f(CUBE[14][0], CUBE[14][1], CUBE[14][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[7][0], CUBE[7][1], CUBE[7][2]);
        glVertex3f(CUBE[15][0], CUBE[15][1], CUBE[15][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[14][0], CUBE[14][1], CUBE[14][2]);
        glVertex3f(CUBE[15][0], CUBE[15][1], CUBE[15][2]);
        glEnd();
    }
    Clean();

    //7-ая грань
    if (FACEisReal(X3, X7, X15, X11) == 1) {
        glBegin(GL_LINES);
        glVertex3f(CUBE[3][0], CUBE[3][1], CUBE[3][2]);
        glVertex3f(CUBE[7][0], CUBE[7][1], CUBE[7][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[3][0], CUBE[3][1], CUBE[3][2]);
        glVertex3f(CUBE[11][0], CUBE[11][1], CUBE[11][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[7][0], CUBE[7][1], CUBE[7][2]);
        glVertex3f(CUBE[15][0], CUBE[15][1], CUBE[15][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[11][0], CUBE[11][1], CUBE[11][2]);
        glVertex3f(CUBE[15][0], CUBE[15][1], CUBE[15][2]);
        glEnd();

    }
    Clean();

    //8-ая грань
    if (FACEisReal(X2, X3, X11, X10) == 1) {
        glBegin(GL_LINES);
        glVertex3f(CUBE[2][0], CUBE[2][1], CUBE[2][2]);
        glVertex3f(CUBE[3][0], CUBE[3][1], CUBE[3][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[2][0], CUBE[2][1], CUBE[2][2]);
        glVertex3f(CUBE[10][0], CUBE[10][1], CUBE[10][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[3][0], CUBE[3][1], CUBE[3][2]);
        glVertex3f(CUBE[11][0], CUBE[11][1], CUBE[11][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[10][0], CUBE[10][1], CUBE[10][2]);
        glVertex3f(CUBE[11][0], CUBE[11][1], CUBE[11][2]);
        glEnd();
    }
    Clean();

    //9-ая грань
    if (FACEisReal2(X0, X1, X5, X3, X4, X2, X6, X7) == 1) {
        glBegin(GL_LINES);
        glVertex3f(CUBE[0][0], CUBE[0][1], CUBE[0][2]);
        glVertex3f(CUBE[1][0], CUBE[1][1], CUBE[1][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[1][0], CUBE[1][1], CUBE[1][2]);
        glVertex3f(CUBE[5][0], CUBE[5][1], CUBE[5][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[5][0], CUBE[5][1], CUBE[5][2]);
        glVertex3f(CUBE[4][0], CUBE[4][1], CUBE[4][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[4][0], CUBE[4][1], CUBE[4][2]);
        glVertex3f(CUBE[6][0], CUBE[6][1], CUBE[6][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[6][0], CUBE[6][1], CUBE[6][2]);
        glVertex3f(CUBE[7][0], CUBE[7][1], CUBE[7][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[7][0], CUBE[7][1], CUBE[7][2]);
        glVertex3f(CUBE[3][0], CUBE[3][1], CUBE[3][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[3][0], CUBE[3][1], CUBE[3][2]);
        glVertex3f(CUBE[2][0], CUBE[2][1], CUBE[2][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[2][0], CUBE[2][1], CUBE[2][2]);
        glVertex3f(CUBE[0][0], CUBE[0][1], CUBE[0][2]);
        glEnd();
    }
    Clean();

    //10-ая грань
    if (FACEisReal2(X8, X9, X13, X11, X12, X10, X14, X15) == 1) {
        glBegin(GL_LINES);
        glVertex3f(CUBE[8][0], CUBE[8][1], CUBE[8][2]);
        glVertex3f(CUBE[9][0], CUBE[9][1], CUBE[9][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[9][0], CUBE[9][1], CUBE[9][2]);
        glVertex3f(CUBE[13][0], CUBE[13][1], CUBE[13][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[13][0], CUBE[13][1], CUBE[13][2]);
        glVertex3f(CUBE[12][0], CUBE[12][1], CUBE[12][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[12][0], CUBE[12][1], CUBE[12][2]);
        glVertex3f(CUBE[14][0], CUBE[14][1], CUBE[14][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[14][0], CUBE[14][1], CUBE[14][2]);
        glVertex3f(CUBE[15][0], CUBE[15][1], CUBE[15][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[15][0], CUBE[15][1], CUBE[15][2]);
        glVertex3f(CUBE[11][0], CUBE[11][1], CUBE[11][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[11][0], CUBE[11][1], CUBE[11][2]);
        glVertex3f(CUBE[10][0], CUBE[10][1], CUBE[10][2]);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(CUBE[10][0], CUBE[10][1], CUBE[10][2]);
        glVertex3f(CUBE[8][0], CUBE[8][1], CUBE[8][2]);
        glEnd();
    }
    Clean();

    glFlush();
    glutSwapBuffers();

}

void resize(int width, int height) {
    glutReshapeWindow(600, 600);
}
// Обработка нажимаемых пользователем клавиш
void specialKeys(int key, int x, int y) {
    switch (key) {
    case GLUT_KEY_RIGHT:
        rotate_y += 5;
        break;
    case GLUT_KEY_LEFT:
        rotate_y -= 5;
        break;
    case GLUT_KEY_UP:
        rotate_x += 5;
        break;
    case GLUT_KEY_DOWN:
        rotate_x -= 5;
        break;
    case GLUT_KEY_PAGE_UP:
        rotate_z += 5;
        break;
    case GLUT_KEY_PAGE_DOWN:
        rotate_z -= 5;
        break;
    case GLUT_KEY_HOME:
        transf = 0.9;
        break;
    case GLUT_KEY_END:
        transf = 1.1;
        break;
    case GLUT_KEY_F1:
        OZ = 1;
        break;
    case GLUT_KEY_F2:
        OZ2 = 1;
        break;
    case GLUT_KEY_F3:
        OX = 1;
        break;
    }

    glutPostRedisplay();

}

int main(int argc, char* argv[]) {
    glutInit(&argc, argv);
    glutInitWindowSize(600, 600);
    glutCreateWindow("Lab2");
    glutDisplayFunc(display);
    glutSpecialFunc(specialKeys);
    glutReshapeFunc(resize);
    glutMainLoop();
    return 0;
}