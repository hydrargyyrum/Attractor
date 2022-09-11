#include <iostream>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <random>
#include <conio.h>
#include "attractor_ui.h"

using namespace std;
const double g = 9.8;
const int m = 2;
const int a = 1;
//const int n = 5;
const int n = 12;
const double Cx = 0.47;
const double S = 0.00785;
const double mass = 1;


class gauss {
private:
	double mean;
	double stddev;

	double uniformRand()
	{
		return static_cast<double>(rand()) / RAND_MAX;
	}

public:

	double normalRand()
	{
		return sqrt(-2.0 * log(uniformRand())) * cos(2.0 * 3.14 * uniformRand()) * stddev + mean;
	}

	gauss(double mean, double stddev) {
		this->mean = mean;
		this->stddev = stddev;
	}
};

class IMathModel {
protected:
	vector<double> X;
	vector<vector<double>> Vectors;
public:
	virtual vector<double> RigthParts() = 0;

	virtual void SetX0(vector<double> x0) {
		X = x0;
		Vectors.push_back(x0);
	}

	virtual double VGX() { return 0; }
	virtual double VGY() { return 0; }

	vector<double> X0() { return X; }
	vector<vector<double>> X_all() { return Vectors; }

	int size() { return X.size(); }

	void print() {
		int n = X.size();
		for (int i = 0; i < n; i++) {
			cout << X[i] << " ";
		}
		cout << endl;
	};
};



class MotionModel : public IMathModel {
protected:
	double VGx;
	double VGy;
	gauss* Gauss;
public:
	MotionModel(vector<double>m, gauss* Gauss) {
		this->Gauss = Gauss;
		SetX0(m);
		VGx = Gauss->normalRand();
		VGy = Gauss->normalRand();
	}

	double VGX() override { return VGx; }
	double VGY() override { return VGy; }

	vector<double> RigthParts() override {
		vector<double> vTmp(n);
		double x = X[0];
		double y = X[1];
		double Vx = X[2];
		double Vy = X[3];
		vTmp[0] = Vx;
		vTmp[1] = Vy;
		vTmp[2] = -sqrt(Vx * Vx + VGx * VGx) * (Vx - VGx) * Cx * S * 0.5 / mass;
		vTmp[3] = -sqrt(Vy * Vy + VGy * VGy) * (Vy - VGy) * Cx * S * 0.5 / mass - g;
		return vTmp;
	}

	void SetX0(vector<double> x0) override {
		X = x0;
	}
};


class attractor : public IMathModel {
public:
	attractor(vector<double> m) {
		SetX0(m);
	}

	//vector<double> RigthParts() override
	//{
	//	//int n = SizeX();
	//	vector<double> vTmp(n);
	//	double X2 = X[0]; //x
	//	double Vx = X[1];//Vx
	//	double Y = X[2]; //y
	//	double Vy = X[3];//Vy
	//	double Z = X[4]; //z
	//	double Vz = X[5];//Vz
	//	double F = X[6]; //fi
	//	double dF = X[7];//dif fi
	//	double T = X[8]; //teta
	//	double dT = X[9];//dif teta
	//	double P = X[10];//psi
	//	double dP = X[11];//dif psi
	//	double X1 = 0, Y1 = 0, Z1 = 0, F1 = 0, T1 = 0, P1 = 0;
	//	double Ix = cos(T) * sin(P);
	//	double Iy = cos(F) * sin(P);
	//	double Iz = cos(F) * cos(T);
	//	double Ir = sin(F) * sin(T) * sin(P);
	//	vTmp[0] = X1;
	//	vTmp[1] = -(cos(F) * sin(T) * cos(P) + sin(F) * sin(P)) * 1000.0 / 10 ;
	//	vTmp[2] = Y1;
	//	vTmp[3] = -(cos(F) * sin(T) * sin(P) + sin(F) * cos(P)) * 1000.0 / 10;
	//	vTmp[4] = Z1;
	//	vTmp[5] = g - (cos(F) * cos(T)) * 1000.0 / 10 ;
	//	vTmp[6] = F1;
	//	vTmp[7] = T1 * P1 * (Iy - Iz) / Ix - Ir / Ix * T1 * 9.81;
	//	vTmp[8] = T1;
	//	vTmp[9] = F1 * P1 * (Iz - Ix) / Iy - Ir / Iy * F1 * 9.81;
	//	vTmp[10] = P1;
	//	vTmp[11] = T1 * F1 * (Ix - Iy) / Iz;
	//	return vTmp;
	//}
	vector<double> RigthParts() override {
		vector<double> vTmp(n);
		double x = X[0];
		double y = X[1];
		double z = X[2];
		double Vy = X[3];

		vTmp[0] = 10 * (y - x);
		vTmp[1] = x * (28 - z) - y;
		vTmp[2] = x * y - 8.0 / 3.0 * z;
		vTmp[3] = 0;
		return vTmp;
	}

	//void SetX0(vector<double> x0) override {
	//	X = x0;
	//}
};

struct integratorArgs {
	integratorArgs(IMathModel* m, double tk) {
		model = m;
		this->tk = tk;
	}
	IMathModel* model = nullptr;
	double tk;
};

struct eylerArgs : integratorArgs {
	eylerArgs(IMathModel* m, double tk, double e, double dt, double t0) :integratorArgs(m, tk), e(e), dt(dt), t0(t0) {};
	double e;
	double dt;
	double t0;
};

struct rungeArgs : integratorArgs {
	rungeArgs(IMathModel* m, double tk, double h) : integratorArgs(m, tk), h(h) {};
	double h;
};


struct eylerCorrectionArgs : integratorArgs {
	eylerCorrectionArgs(IMathModel* m, double tk, double h, double dt, double t0) :integratorArgs(m, tk), h(h), dt(dt), t0(t0) {};
	double dt;
	double h;
	double t0;
};

class integrator {
protected:
	IMathModel* model;
	double tk;
public:

	integrator(integratorArgs args) : tk(args.tk), model(args.model) { }

	~integrator() {
		cout << "class is destroyed" << endl;
	}

	vector<double> RigthParts() {
		return model->RigthParts();
	}

	void SetX0(vector<double> x0) {
		model->SetX0(x0);
		//cout << "X = " << x0[0] << " Y = " << x0[1] << " Vx = " << x0[2] << endl;
	}

	void print() {
		model->print();
	}

	virtual void integration() = 0;
};

class eyler : public integrator {
protected:
	ofstream* fout;
	const double e;
	const double dt;
	const double tk;
	const double t0;
public:

	void printF(vector<double> x0, double t) {
		ofstream& pFout = *fout;
		pFout << "   t=" << t << "    x=" << x0[0] << "    y=" << x0[1] << endl;
	}

	eyler(eylerArgs args) : integrator(args), dt(args.dt), e(args.e), tk(args.tk), t0(args.t0) {
		fout = new ofstream("f.txt");
	}

	~eyler() {
		fout->close();
		delete fout;
	}

	double t0_() { return t0; }

	void integration() override {
		double t = t0_();
		int n = model->size();
		vector<double> vPrev = model->X0();
		double y = vPrev[1];
		double x = vPrev[0];
		while (t < tk) {
			vector<double> f = RigthParts();
			vector<double> vNext(n);
			t = t + dt;
			for (int i = 0; i < n; i++)
				vNext[i] = vPrev[i] + f[i] * dt;

			/// ATTRACTOR
			SetX0(vNext);
			printF(vNext, t);
			vPrev = vNext;

			/// VECTOR
			 /* if (vNext[1] >= 0)
			{
				SetX0(vNext);
				printF(vNext, t);
				vPrev = vNext;
				y = vNext[1];
			}
			else
				break; */
		}
		//cout << "S=" << model->X0()[0] - x << endl; 
	}
};

class runge : public integrator {
protected:
	ofstream* fout;
	const double h;
public:
	runge(rungeArgs args) : integrator(args), h(args.h) {
		fout = new ofstream("f.txt");
	}

	~runge() {
		fout->close();
		delete fout;
	}

	void printF(vector<double> x0, double t) {
		ofstream& pFout = *fout;
		pFout << "   t=" << t << "    x=" << x0[0] << "    y=" << x0[1] << endl;
	}

	struct Point {
		double x, y, vx, vy;

		Point() {}
		Point(double x, double y, double vx, double vy) : x(x), y(y), vx(vx), vy(vy) { }
		Point add(double d) {
			return Point(x + d, y + d, vx + d, vy + d);
		}
	};

	double f(int i, double t, Point point) {
		switch (i) {
		case 0:
			/// ATTRACTOR
			return 10 * (point.y - point.x);

			/// VECTOR
			//return point.vx;
		case 1:
			/// ATTRACTOR
			return point.x * (28 - point.vx) - point.y;

			/// VECTOR
			//return point.vy;
		case 2: {
			/// ATTRACTOR
			return point.x * point.y - 8.0 / 3.0 * point.vx;

			/// VECTOR
			//double Vx = point.vx;
			//double VGx = this->model->VGX();
			//return -sqrt(Vx * Vx + VGx * VGx) * (Vx - VGx) * Cx * S * 0.5 / mass;
		}
		case 3:
			double Vy = point.vy;
			double VGy = this->model->VGY();
			return -sqrt(Vy * Vy + VGy * VGy) * (Vy - VGy) * Cx * S * 0.5 / mass - g;
		}
	}

	double oneStep(int i, vector<double> vPrev, double t) {
		Point p;
		p.x = vPrev[0];
		p.y = vPrev[1];
		p.vx = vPrev[2];
		p.vy = vPrev[3];

		double k1 = h * f(i, t, p);
		double k2 = h * f(i, t, p.add(k1 / 2));
		double k3 = h * f(i, t, p.add(k2 / 2));
		double k4 = h * f(i, t, p.add(k3));
		double res = (vPrev[i] + 1.0 / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4));
		return res;
	}

	void integration() override {
		int n_count = 1000;
		int n = model->size();
		double x = model->X0()[0];
		double t = 0;
		vector<double> vPrev = model->X0();
		for (int i = 0; i < n_count; i++) {
			t = t + 1;
			vector<double> vNext(n);
			for (int j = 0; j < n; j++) {
				vNext[j] = oneStep(j, vPrev, t);
			}
			vPrev = vNext;
			/// ATTRACTOR
			SetX0(vNext);
			printF(vNext, t);
			vPrev = vNext;

			/// VECTOR
			 /* if (vNext[1] >= 0)
			{
				SetX0(vNext);
				printF(vNext, t);
				vPrev = vNext;
				y = vNext[1];
			}
			else
				break; */
		}
	}
};



class eylerCorrection : public integrator {
protected:
	ofstream* fout;
	const double h;
	const double dt;
	const double tk;
	const double t0;
public:

	eylerCorrection(eylerCorrectionArgs args) : integrator(args), dt(args.dt), h(args.h), tk(args.tk), t0(args.t0) { fout = new ofstream("f.txt"); }
	~eylerCorrection() {
		fout->close();
		delete fout;
	}

	void printF(vector<double> x0, double t) {
		ofstream& pFout = *fout;
		pFout << "   t=" << t << "    x=" << x0[0] << "    y=" << x0[1] << endl;
	}

	double t0_() { return t0; }

	struct Point {
		double x, y, vx, vy;

		Point() {}
		Point(double x, double y, double vx, double vy) : x(x), y(y), vx(vx), vy(vy) { }

		Point addX(double d) {
			return Point(x + d, y, vx, vy);
		}
		Point addY(double d) {
			return Point(x, y + d, vx, vy);
		}
		Point addVx(double d) {
			return Point(x, y, vx + d, vy);
		}
		Point addVy(double d) {
			return Point(x, y, vx, vy + d);
		}
	};



	double f(int i, Point point) {
		switch (i) {
		case 0:
			return point.vx;
		case 1:
			return point.vy;
		case 2: {
			double Vx = point.vx;
			double VGx = this->model->VGX();
			return -sqrt(Vx * Vx + VGx * VGx) * (Vx - VGx) * Cx * S * 0.5 / mass;
		}
		case 3:
			double Vy = point.vy;
			double VGy = this->model->VGY();
			return -sqrt(Vy * Vy + VGy * VGy) * (Vy - VGy) * Cx * S * 0.5 / mass - g;
		}
	}


	double oneStep(int i, vector<double> vPrev) {
		Point p;
		p.x = vPrev[0];
		p.y = vPrev[1];
		p.vx = vPrev[2];
		p.vy = vPrev[3];
		switch (i) {
		case 0:
			return (p.x + h / 2 * (f(i, p) + f(i, p.addX(h * f(i, p)))));
		case 1:
			return (p.y + h / 2 * (f(i, p) + f(i, p.addY(h * f(i, p)))));
		case 2:
			return (p.vx + h / 2 * (f(i, p) + f(i, p.addVx(h * f(i, p)))));
		case 3:
			return (p.vy + h / 2 * (f(i, p) + f(i, p.addVy(h * f(i, p)))));
		}
	}


	void integration() override {
		int n = model->size();
		double t = t0_();
		double x = model->X0()[0];
		vector<double> vPrev = model->X0();
		while (t < tk) {
			t = t + dt;
			vector<double> vNext(n);
			for (int j = 0; j < n; j++) {
				vNext[j] = oneStep(j, vPrev);
			}
			if (vNext[1] < 0)
				break;
			else
			{
				SetX0(vNext);
				printF(vNext, t);
				vPrev = vNext;
			}
		}
		cout << "S=" << (model->X0()[0] - x) << endl;

	}
};


// ctrl+shift+b - build project
// ctrl+k+c - comment
// ctrl+k+u - uncomment
// ctrl+k+d - auto form
// BL don't combine with UI!!!!!!!!

void main() {
	//{
	//	attractor_show_args args_attr;
	//	args_attr.example_to_show = true;
	//	args_attr.example_continiously_draw = false;
	//	attractor_show(args_attr);
	//	return;
	//}

	srand(time(nullptr));
	gauss G(-1, 4);
	vector<double> m = { 10, 10, 10, 3 };
	//vector<double> m = { 10, 10, 10,10,10,10,10,10,10,10,10,10};
	//MotionModel* model = new MotionModel(m, &G);

	attractor* modelA = new attractor(m);
	eylerArgs args1(modelA, 10, 0.01, 0.01, 0);
	eyler c1(args1);

	/*eylerCorrectionArgs args2(model, 100, 0.01, 1, 0);
	eylerCorrection c1(args2); */

	/*rungeArgs args(modelA, 1000, 0.01);
	runge c1(args);*/

	c1.integration();
	//cin.get();

	{
		//return;
		attractor_show_args args_attr;
		args_attr.example_to_show = false;
		args_attr.example_continiously_draw = true;
		args_attr.points = modelA->X_all();
		attractor_show(args_attr);
	}
}