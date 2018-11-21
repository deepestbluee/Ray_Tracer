#ifndef __HW1__HELPER__
#define __HW1__HELPER__

#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <math.h>
#include <limits>

float scalarVector(const parser::Vec3f &lhs ,const parser::Vec3f &rhs);

int scalarVector(const parser::Vec3i &lhs ,const parser::Vec3i &rhs);

parser::Vec3f crossVector(const parser::Vec3f &lhs ,const parser::Vec3f &rhs);

parser::Vec3i crossVector(const parser::Vec3i &lhs ,const parser::Vec3i &rhs);

parser::Vec3i sum_vector(const parser::Vec3i &a,const parser::Vec3i &b);

parser::Vec3f sum_vector(const parser::Vec3f &a,const parser::Vec3f &b);

parser::Vec3i subtract_vector(const parser::Vec3i &a,const parser::Vec3i &b);

parser::Vec3f subtract_vector(const parser::Vec3f &a,const parser::Vec3f &b);

float normVector(const parser::Vec3f &a);

float distanceVector(const parser::Vec3f &a, const parser::Vec3f &b);

parser::Vec3f divideVector(const float &lhs,const parser::Vec3f &rhs);


template<class T>
parser::Vec3f multiVector(const T &lhs,const parser::Vec3f &rhs)
{
    parser::Vec3f c;
    c.x = lhs*rhs.x;
    c.y = lhs*rhs.y;
    c.z = lhs*rhs.z;
    return c;  
}

template<class T>
parser::Vec3i multiVector(const T &lhs,const parser::Vec3i &rhs)
{
    parser::Vec3i c;
    c.x = lhs*rhs.x;
    c.y = lhs*rhs.y;
    c.z = lhs*rhs.z;
    return c;  
}

parser::Vec3f computeRayEquation(const int &y,
                                 const int &x,
                                 const parser::Camera &cam);

float spheresIntersect(const parser::Vec3f &ray,
                       const parser::Vec3f &pos,
                       const parser::Vec3f &sphereCenter,
                       const float &rad);

float determinant(const float &a,
                  const float &b,
                  const float &c,
                  const float &d,
                  const float &e,
                  const float &f,
                  const float &g,
                  const float &h,
                  const float &i);

float trianglesIntersect(const parser::Vec3f &ray,
                          const parser::Vec3f &position,
                          const std::vector<parser::Vec3f> &vertex_data,
                          const parser::Triangle &triangle);

float meshHelper(const parser::Vec3f &ray,
                          const parser::Vec3f &position,
                          const std::vector<parser::Vec3f> &vertex_data,
                          const parser::Face &face);

std::vector<float> meshIntersect(const parser::Vec3f &ray,
                          const parser::Vec3f &position,
                          const std::vector<parser::Vec3f> &vertex_data,
                          const parser::Mesh &mesh);

parser::Vec3f returnS(const int &y,
                                 const int &x,
                                 const parser::Camera &cam);

std::vector<int> mirrorRecur(const parser::Scene &scene,
                          const parser::Vec3f &cam_pos,
                          const parser::Vec3f &ray,
                          const int depth);
#endif
