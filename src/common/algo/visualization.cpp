// #include "common/algo/visualization.hpp"
//
// #include <iostream>
// #include <sstream>
// #include <string>
// #include <vector>
//
// #ifdef _WIN32
// #else
// #include <sys/stat.h>  // For mkdir in POSIX
// #include <sys/time.h>  // For gettimeofday
// #include <sys/types.h> // For types in mkdir
//
// #include <Eigen/Core>
// #include <matplot/matplot.h>
//
// using namespace std;
// using namespace matplot;
// using namespace Eigen;
//
// #define LOCATION "/home/"
//
// namespace common {
// namespace algo {
// std::string visualize_icp_iteration(const MatrixX2d &_X, const MatrixX2d &Y,
//                                     const MatrixX2d &X,
//                                     const MatrixX2d &samples,
//                                     const MatrixX2d &C, int iteration,
//                                     const MatrixXd &R, const VectorXd &t) {
//   if (!VISUALIZE) {
//     return "";
//   }
//   struct timeval tv;
//   gettimeofday(&tv, NULL);
//   struct tm *ltm = localtime(&tv.tv_sec);
//
//   const char *home = std::getenv("HOME");
//   std::string homeDir = home ? std::string(home) : "";
//
//   std::stringstream newDirName;
//   newDirName << homeDir << LOCATION << 1900 + ltm->tm_year << std::setw(2)
//              << std::setfill('0') << 1 + ltm->tm_mon << std::setw(2)
//              << std::setfill('0') << ltm->tm_mday << '-' << std::setw(2)
//              << std::setfill('0') << ltm->tm_hour << std::setw(2)
//              << std::setfill('0') << ltm->tm_min << std::setw(2)
//              << std::setfill('0') << ltm->tm_sec << std::setw(3)
//              << std::setfill('0') << tv.tv_usec / 1000;
//   const std::string newDirNameStr = newDirName.str();
//
//   if (mkdir(newDirNameStr.c_str(), 0777) == -1) {
//     std::cerr << "Error creating directory '" << newDirNameStr
//               << "': " << strerror(errno) << std::endl;
//   }
//
//   visualize_icp_iteration(_X, Y, X, samples, C, iteration, R, t,
//   newDirNameStr); return newDirNameStr;
// }
//
// /**
//  * @brief Visualizes the query points, target points, transformed points.
//  *
//  * @param _X
//  * @param Y
//  * @param X
//  * @param samples
//  * @param C
//  * @param iteration
//  */
// void visualize_icp_iteration(const MatrixX2d &_X, const MatrixX2d &Y,
//                              const MatrixX2d &X, const MatrixX2d &samples,
//                              const MatrixX2d &C, int iteration,
//                              const MatrixXd &R, const VectorXd &t,
//                              std::string dirname) {
//   if (!VISUALIZE) {
//     return;
//   }
//
//   // Create a new figure
//   auto fig = figure(true);
//   fig->size(1200, 1200); // Set the size of the figure (in pixels
//   auto ax = fig->current_axes();
//   // ax->axes_aspect_ratio(1.0);
//
//   hold(ax, true);
//   // Plot original points _X in red
//   vector<double> x1(_X.col(0).data(), _X.col(0).data() + _X.rows());
//   vector<double> y1(_X.col(1).data(), _X.col(1).data() + _X.rows());
//   ax->plot(x1, y1, "r-o")
//       ->marker_size(4)
//       .marker_face_color("red")
//       .marker_face_alpha(0.9);
//
//   // Plot target points Y in blue
//   vector<double> x2(Y.col(0).data(), Y.col(0).data() + Y.rows());
//   vector<double> y2(Y.col(1).data(), Y.col(1).data() + Y.rows());
//   ax->plot(x2, y2, "b-o")
//       ->marker_size(4)
//       .marker_face_color("blue")
//       .marker_face_alpha(0.9);
//
//   // Plot transformed points X in green
//   vector<double> x3(X.col(0).data(), X.col(0).data() + X.rows());
//   vector<double> y3(X.col(1).data(), X.col(1).data() + X.rows());
//   ax->plot(x3, y3, "g-o")
//       ->marker_size(4)
//       .marker_face_color("green")
//       .marker_face_alpha(0.9);
//
//   // Plot sample points and corresponding points in matching colors with
//   // lines
//   vector<double> x4(samples.col(0).data(),
//                     samples.col(0).data() + samples.rows());
//   vector<double> y4(samples.col(1).data(),
//                     samples.col(1).data() + samples.rows());
//   vector<double> x5(C.col(0).data(), C.col(0).data() + C.rows());
//   vector<double> y5(C.col(1).data(), C.col(1).data() + C.rows());
//
//   for (long i = 0; i < samples.rows(); ++i) {
//     ax->scatter(vector<double>{x4[i]},
//     vector<double>{y4[i]})->marker_size(5);
//     ax->scatter(vector<double>{x5[i]},
//     vector<double>{y5[i]})->marker_size(5);
//
//     ax->plot({x4[i], x5[i]}, {y4[i], y5[i]})
//         ->line_width(1)
//         .marker_size(1)
//         .color({0.5, 0.5, 0.5, 0.5});
//   }
//
//   // Display R matrix and t vector
//   // Display in the top left of the plot
//   auto x_lim = ax->xlim();
//   auto y_lim = ax->ylim();
//
//   stringstream rt_info;
//   double theta = acos(R(0, 0));
//   rt_info << "R\\n"
//           << "theta: " << theta << "\\n";
//
//   for (int i = 0; i < R.rows(); ++i) {
//     for (int j = 0; j < R.cols(); ++j) {
//       rt_info << R(i, j) << " ";
//     }
//     rt_info << "\\n";
//   }
//
//   rt_info << "t:\\n";
//   for (int i = 0; i < t.size(); ++i) {
//     rt_info << t(i) << "\\n";
//   }
//   ax->text(x_lim[0], y_lim[1], rt_info.str());
//   //   ->alignment(labels::alignment::left);
//
//   hold(ax, false);
//
//   // Set labels and title
//   ax->axis(matplot::equal);
//   ax->xlabel("X-axis");
//   ax->ylabel("Y-axis");
//
//   stringstream title;
//   title << "ICP Iteration: " << iteration << " Samples: " << samples.rows();
//   ax->title(title.str());
//
//   //   Save the plot as an image
//   stringstream filename;
//   filename << dirname << "_icp_iteration_" << setw(4) << setfill('0')
//            << iteration << ".png";
//   fig->save(filename.str());
// }
//
// } // namespace algo
// } // namespace common
// #endif
