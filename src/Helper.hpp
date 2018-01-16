#ifndef HELPER_HPP
#define HELPER_HPP

#include <math.h>
#include <iostream>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"

using namespace std;

class Helper{
  public:

    // For converting back and forth between radians and degrees.
    static double deg2rad(double x) { return x * M_PI / 180; }
    static double rad2deg(double x) { return x * 180 / M_PI; }
    static double kmh2mps(double x) { return x / 3.6; }
    static double mph2mps(double x) { return x / 2.24; }

    // Checks if the SocketIO event has JSON data.
    // If there is data the JSON object in string format will be returned,
    // else the empty string "" will be returned.
    static string hasData(string s) {
      auto found_null = s.find("null");
      auto b1 = s.find_first_of("[");
      auto b2 = s.rfind("}]");
      if (found_null != string::npos) {
        return "";
      } else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
      }
      return "";
    }

    // Evaluate a polynomial.
    static double polyeval(Eigen::VectorXd coeffs, double x) {
      double result = 0.0;
      for (int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, i);
      }
      return result;
    }

    static vector<Eigen::VectorXd> convertCoordinates(vector<double> ptsx, vector<double> ptsy, double px, double py, double psi){
      Eigen::VectorXd state(6);
      int waypoints_size = ptsx.size();
      Eigen::VectorXd waypoints_x(waypoints_size);
      Eigen::VectorXd waypoints_y(waypoints_size);

      for(int i=0; i<waypoints_size; ++i){
        double diff_x = (ptsx[i] - px);
        double diff_y = (ptsy[i] - py);
        waypoints_x[i] = (diff_x * cos(-psi)) - (diff_y * sin(-psi));
        waypoints_y[i] = -(diff_x * sin(-psi)) - (diff_y * cos(-psi));
      }

      return {waypoints_x, waypoints_y};
    }

    // Fit a polynomial.
    // Adapted from
    // https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
    static Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                            int order) {
      assert(xvals.size() == yvals.size());
      assert(order >= 1 && order <= xvals.size() - 1);
      Eigen::MatrixXd A(xvals.size(), order + 1);

      for (int i = 0; i < xvals.size(); i++) {
        A(i, 0) = 1.0;
      }

      for (int j = 0; j < xvals.size(); j++) {
        for (int i = 0; i < order; i++) {
          A(j, i + 1) = A(j, i) * xvals(j);
        }
      }

      auto Q = A.householderQr();
      auto result = Q.solve(yvals);
      return result;
    }
};

#endif /* HELPER_HPP */