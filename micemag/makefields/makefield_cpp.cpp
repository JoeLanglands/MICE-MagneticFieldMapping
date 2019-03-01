#include <boost/python.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <cmath>


#define MU_0 M_PI*4e-7


/* Functional code because making it into a class that is interfaceable with python is a massive
   sweffort.  It will instead be used and wrapped into a python class for ease of use.
 */

double ellipk(double _arg){
  return boost::math::ellint_1(_arg);
}

double ellipe(double _arg){
  return boost::math::ellint_2(_arg);
}


std::tuple<double, double> field_from_loop(double current, double radius, double R, double Z){
  double rho = sqrt(R*R + Z*Z);
  double alpha = sqrt(radius*radius + rho*rho - 2*radius*R);
  double beta = sqrt(radius*radius + rho*rho + 2*radius*R);
  double k = sqrt(1 - (alpha*alpha)/(beta*beta));

  double C = (MU_0*current)/M_PI;
  double a = 1.0/(2*(alpha*alpha)*beta);
  double b = Z/(2*(alpha*alpha)*beta*R);

  double Bz = C*a*((radius*radius - rho*rho)*ellipe(k) + (alpha*alpha)*ellipk(k));
  double Br = C*b*((radius*radius + rho*rho)*ellipe(k) - (alpha*alpha)*ellipk(k));

  if (std::isnan(Bz) == true || std::isinf(Bz) == true){
    Bz = 0;
  }
  if (std::isnan(Br) == true || std::isinf(Br) == true){
    Br = 0;
  }
  
  return std::make_tuple(Br, Bz);
}



std::tuple<double, double> make_current_layer(int nTurns, double separation, double startPos, double current, double radius, double R, double Z){

  double offset = startPos - separation;
  double newZ = 0;
  std::tuple<double, double> newBrBz (0, 0);
  double totalBr = 0;
  double totalBz = 0;

  for (int i = 0; i < nTurns; i++){
    offset += separation;
    newZ = Z - offset;
    newBrBz = field_from_loop(current, radius, R, newZ);
    totalBr += std::get<0>(newBrBz);
    totalBz += std::get<1>(newBrBz);
  }
  return std::make_tuple(totalBr, totalBz);

}


std::tuple<double, double> make_coil_from_layers(int nLayers, int nTurns, double layerSep, double loopSep, double startPos, double current, double minR, double R, double Z){
  double totalBr = 0;
  double totalBz = 0;

  std::tuple<double, double> newBrBz (0, 0);

  double radius = minR - layerSep;

  for (int i = 0; i < nLayers; i++){
    radius += layerSep;
    newBrBz = make_current_layer(nTurns, loopSep, startPos, current, radius, R, Z);
    totalBr += std::get<0>(newBrBz);
    totalBz += std::get<1>(newBrBz);
  }
  return std::make_tuple(totalBr, totalBz);
}


std::tuple<double, double> make_coil(double currentDensity, double centrePos, double length, double rInner, double rOuter, int nLayers, int nTurns, double R, double Z){
  double z1 = centrePos - 0.5*length;
  double z2 = centrePos + 0.5*length;
  double height = rOuter - rInner;

  double loopSep = length/double(nTurns);
  double layerSep = height/double(nLayers);

  double firstLoopZ = z1 + 0.5*loopSep;
  double firstLoopR = rInner + 0.5*layerSep;

  int totalLoops = nLayers*nTurns;
  double totalCurrent = length*height*currentDensity*1.0e6; //convert from A/mm^2 to A/m^2
  double currentPerLoop = totalCurrent/double(totalLoops);

  std::tuple<double, double> newBrBz = make_coil_from_layers(nLayers, nTurns, layerSep, loopSep, firstLoopZ, currentPerLoop, firstLoopR, R, Z);

  return newBrBz;
}

double current_to_density(double current, double rInner, double rOuter, double length, int nLayers, int nTurns){
  return (current*double(nLayers)*double(nTurns))/(length*(rOuter - rInner)*1e6); //in A/mm^2
}

std::tuple<double, double, double, double, double, double> convert_to_cartesian(double r, double phi, double z, double Br, double Bphi, double Bz){
  using namespace boost::numeric::ublas;
  double x = r*cos(phi);
  double y = r*sin(phi);
  
  //form vector of (Br, Bphi)
  vector<double> B(2);
  B(0) = Br; B(1) = Bphi;
  
  //form matrix that converts B_r/phi to B_x/y
  matrix<double> xy_rphi(2, 2);
  xy_rphi(0,0) = cos(phi); xy_rphi(0,1) = -1.0*sin(phi);
  xy_rphi(1,0) = sin(phi); xy_rphi(1,1) = cos(phi);
  
  vector<double> Bxy = prod(xy_rphi, B);
  
  return std::make_tuple(x, y, z, Bxy(0), Bxy(1), Bz);
}


std::tuple<double, double, double, double, double, double> convert_to_polar(double x, double y, double z, double Bx, double By, double Bz){
  using namespace boost::numeric::ublas;
  double r = sqrt(x*x + y*y);
  double phi = atan2(y, x);

  //form vector of (Bx, By)
  vector<double> B(2);
  B(0) = Bx; B(1) = By;

  //form matrix that converts B_x/y to B_r/phi
  matrix<double> rphi_xy(2, 2);
  rphi_xy(0, 0) = cos(phi); rphi_xy(0, 1) = sin(phi);
  rphi_xy(1, 0) = -1.0*sin(phi) ; rphi_xy(1, 1) = cos(phi);

  vector<double> Brphi = prod(rphi_xy, B);

  return std::make_tuple(r, phi, z, Brphi(0), Brphi(1), Bz);
}

std::tuple<double, double, double, double, double, double> rotate(double x, double y, double z, double Bx, double By, double Bz, double thetaX, double thetaY, double px, double py, bool transpose){
  using namespace boost::numeric::ublas;

  //set up vectors for x,y,z and Bxyz to manipulate on
  vector<double> X(3);
  X(0) = x; X(1) = y; X(2) = z;

  vector<double> B(3);
  B(0) = Bx; B(1) = By; B(2) = Bz;

  //set up offset vector
  vector<double> offset(3);
  offset(0) = px; offset(1) = py; offset(2) = 0;
  
  //set up rotation matrix
  matrix<double> R(3, 3);
  R(0, 0) = cos(thetaY); R(0, 1) = 0; R(0, 2) = sin(thetaY);
  R(1, 0) = sin(thetaX)*sin(thetaY); R(1, 1) = cos(thetaX); R(1, 2) = -1.0*cos(thetaY)*sin(thetaX);
  R(2, 0) = -1.0*sin(thetaY)*cos(thetaX); R(2, 1) = sin(thetaX); R(2, 2) = cos(thetaY)*cos(thetaX);

  if (transpose == false){
    vector<double> newX = prod(R, X) + offset;
    vector<double> newB = prod(R, B);
    return std::make_tuple(newX(0), newX(1), newX(2), newB(0), newB(1), newB(2));
    
  } else{
    vector<double> newX = prod(trans(R), X) - offset;
    vector<double> newB = prod(trans(R), B);
    return std::make_tuple(newX(0), newX(1), newX(2), newB(0), newB(1), newB(2));
  }
}


boost::python::tuple get_field_at_point(double current, double centrePos, double length, double rInner, double rOuter, int nLayers, int nTurns, double thetaX, double thetaY, double px, double py, double R, double PHI, double Z){
  double currentDensity = current_to_density(current, rInner, rOuter, length, nLayers, nTurns);

  //convert coords to cartesian
  auto tmp = convert_to_cartesian(R, PHI, Z, 0, 0, 0);
  
  //rotate with R
  tmp = rotate(std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp), 0, 0, 0, thetaX, thetaY, px, py, false);
  
  //convert back to polar
  tmp = convert_to_polar(std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp), std::get<3>(tmp), std::get<4>(tmp), std::get<5>(tmp));
  
  //now calculate the field
  std::tuple<double, double> newBrBz = make_coil(currentDensity, centrePos, length, rInner, rOuter, nLayers, nTurns, std::get<0>(tmp), std::get<2>(tmp));
  
  //convert to cartesian again INCLUDING B FIELD!!
  auto tmp_2 = convert_to_cartesian(std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp), std::get<0>(newBrBz), 0, std::get<1>(newBrBz));
  
  //Rotate with R_t
  tmp_2 = rotate(std::get<0>(tmp_2), std::get<1>(tmp_2), std::get<2>(tmp_2), std::get<3>(tmp_2), std::get<4>(tmp_2), std::get<5>(tmp_2), thetaX, thetaY, px, py, true);
  
  //convert back to polars
  tmp_2 = convert_to_polar(std::get<0>(tmp_2), std::get<1>(tmp_2), std::get<2>(tmp_2), std::get<3>(tmp_2), std::get<4>(tmp_2), std::get<5>(tmp_2));

  //Should check whether the coords are the same as R, PHI, Z
  
  return  boost::python::make_tuple(std::get<3>(tmp_2), std::get<4>(tmp_2), std::get<5>(tmp_2));
}

//This func is like above but takes (x,y,z) as args and returns Bx, By, Bz for g4blgrid field
boost::python::tuple get_field_at_point_xyz(double current, double centrePos, double length, double rInner, double rOuter, int nLayers, int nTurns, double thetaX, double thetaY, double px, double py, double X, double Y, double Z){
  double currentDensity = current_to_density(current, rInner, rOuter, length, nLayers, nTurns);

  auto tmp = rotate(X, Y, Z, 0, 0, 0, thetaX, thetaY, px, py, false);

  tmp = convert_to_polar(std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp), std::get<3>(tmp), std::get<4>(tmp), std::get<5>(tmp));

  std::tuple<double, double> newBrBz =  make_coil(currentDensity, centrePos, length, rInner, rOuter, nLayers, nTurns, std::get<0>(tmp), std::get<2>(tmp));

  auto tmp_2 = convert_to_cartesian(std::get<0>(tmp), std::get<1>(tmp), std::get<2>(tmp), std::get<0>(newBrBz), 0, std::get<1>(newBrBz));

  tmp_2 = rotate(std::get<0>(tmp_2), std::get<1>(tmp_2), std::get<2>(tmp_2), std::get<3>(tmp_2), std::get<4>(tmp_2), std::get<5>(tmp_2), thetaX, thetaY, px, py, true);

  return boost::python::make_tuple(std::get<3>(tmp_2), std::get<4>(tmp_2), std::get<5>(tmp_2));

}


BOOST_PYTHON_MODULE(makefield_cpp){
  using namespace boost::python;
  def ("get_field_at_point", get_field_at_point);
  def ("get_field_at_point_xyz", get_field_at_point_xyz);
}
