#include "caboodle.h"

/* Functions to calculate metric quantities for an orthogonal coordinate system
 * Functions to calculate the scale factors 
 * Functions to calculate the Christoffel symbols
 * Functions to calculate the cell volumes, surface areas, and edge lengths
 */
double scale_h1(double x, double y, double z) {
    return 1.0;
}
double scale_h2(double x, double y, double z) {
    return 1.0;
}
double scale_h3(double x, double y, double z) {
    return 1.0;
}
double scale_dh11(double x, double y, double z) {
    return 0.0;
}
double scale_dh12(double x, double y, double z) {
    return 0.0;
}
double scale_dh13(double x, double y, double z) {
    return 0.0;
}
double scale_dh21(double x, double y, double z) {
    return 0.0;
}
double scale_dh22(double x, double y, double z) {
    return 0.0;
}
double scale_dh23(double x, double y, double z) {
    return 0.0;
}
double scale_dh31(double x, double y, double z) {
    return 0.0;
}
double scale_dh32(double x, double y, double z) {
    return 0.0;
}
double scale_dh33(double x, double y, double z) {
    return 0.0;
}
double scale_dLx(double x, double y, double z, double dx, double dy, double dz) {
    return scale_h1(x,y,z) * dx;
}   
double scale_dLy(double x, double y, double z, double dx, double dy, double dz) {
    return scale_h2(x,y,z) * dy;
}   
double scale_dLz(double x, double y, double z, double dx, double dy, double dz) {
    return scale_h3(x,y,z) * dz;
}   
double scale_Sxy(double x, double y, double z, double dx, double dy, double dz) {
    return scale_h1(x,y,z) * scale_h2(x,y,z) * dx * dy;
}
double scale_Sxz(double x, double y, double z, double dx, double dy, double dz) {
    return scale_h1(x,y,z) * scale_h3(x,y,z) * dx * dz;
}
double scale_Syz(double x, double y, double z, double dx, double dy, double dz) {
    return scale_h2(x,y,z) * scale_h3(x,y,z) * dy * dz;
}

double scale_vol(double x, double y, double z, double dx, double dy, double dz) {
    return scale_dLx(x,y,z,dx,dy,dz) * scale_dLy(x,y,z,dx,dy,dz) * scale_dLz(x,y,z,dx,dy,dz);
}
