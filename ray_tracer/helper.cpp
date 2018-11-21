#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <limits>
#include "helper.h"

float scalarVector(const parser::Vec3f &lhs ,const parser::Vec3f &rhs)
{
    return lhs.x*rhs.x+lhs.y*rhs.y+lhs.z*rhs.z;
}

int scalarVector(const parser::Vec3i &lhs ,const parser::Vec3i &rhs)
{
    return lhs.x*rhs.x+lhs.y*rhs.y+lhs.z*rhs.z;
}

parser::Vec3f crossVector(const parser::Vec3f &lhs ,const parser::Vec3f &rhs)
{
    parser::Vec3f result;   

    result.x = (lhs.y*rhs.z) - (lhs.z*rhs.y);
    result.y = -((lhs.x*rhs.z) - (lhs.z*rhs.x));
    result.z = (lhs.x*rhs.y) - (lhs.y*rhs.x);

    return result; 
}

parser::Vec3i crossVector(const parser::Vec3i &lhs ,const parser::Vec3i &rhs)
{
    parser::Vec3i result;   

    result.x = (lhs.y*rhs.z) - (lhs.z*rhs.y);
    result.y = -((lhs.x*rhs.z) - (lhs.z*rhs.x));
    result.z = (lhs.x*rhs.y) - (lhs.y*rhs.x);

    return result; 

}

parser::Vec3i sum_vector(const parser::Vec3i &a,const parser::Vec3i &b)
{
    parser::Vec3i c;
    c.x = a.x+b.x;
    c.y = a.y+b.y;
    c.z = a.z+b.z;
    return c;
}

parser::Vec3f sum_vector(const parser::Vec3f &a,const parser::Vec3f &b)
{
    parser::Vec3f c;
    c.x = a.x+b.x;
    c.y = a.y+b.y;
    c.z = a.z+b.z;
    return c;
}

parser::Vec3i subtract_vector(const parser::Vec3i &a,const parser::Vec3i &b)
{
    parser::Vec3i c;
    c.x = a.x-b.x;
    c.y = a.y-b.y;
    c.z = a.z-b.z;
    return c;
}

parser::Vec3f subtract_vector(const parser::Vec3f &a,const parser::Vec3f &b)
{
    parser::Vec3f c;
    c.x = a.x-b.x;
    c.y = a.y-b.y;
    c.z = a.z-b.z;
    return c;
}

float normVector(const parser::Vec3f &a)
{
    float c = sqrtf(a.x*a.x+a.y*a.y+a.z*a.z);
    return c;
}

float distanceVector(const parser::Vec3f &a, const parser::Vec3f &b)
{
    float c = sqrtf(powf(a.x-b.x,2)+powf(a.y-b.y,2)+powf(a.z-b.z,2));
    return c;
}

parser::Vec3f divideVector(const float &lhs,const parser::Vec3f &rhs)
{
    parser::Vec3f c;
    c.x = rhs.x/lhs;
    c.y = rhs.y/lhs;
    c.z = rhs.z/lhs;
    return c;  
}

parser::Vec3f computeRayEquation(const int &y, const int &x, const parser::Camera &cam)
{
    float l = cam.near_plane.x;
    float r = cam.near_plane.y;
    float t = cam.near_plane.w;
    float b = cam.near_plane.z;

    parser::Vec3f w = divideVector(normVector(cam.gaze),multiVector(-1,cam.gaze));
    parser::Vec3f v = cam.up;
    parser::Vec3f u = crossVector(v,w);
    u = divideVector(normVector(u),u);
    v = crossVector(w,u);  
    parser::Vec3f m = sum_vector(cam.position , multiVector(cam.near_distance,divideVector(normVector(cam.gaze),cam.gaze)));
    parser::Vec3f q = sum_vector(sum_vector(m,multiVector(l,u)),multiVector(t,v));

    float su = ((r-l)*( x+0.5))/cam.image_width;
    float sv = ((t-b)*( y+0.5))/cam.image_height;

    parser::Vec3f temp1 = multiVector(su,u);
    parser::Vec3f temp2 = multiVector(sv*(-1),v);
    
    parser::Vec3f s = sum_vector(sum_vector(q,temp1),temp2);
    parser::Vec3f result = subtract_vector(s,cam.position);
    
    return result;
}

parser::Vec3f returnS(const int &y, const int &x, const parser::Camera &cam)
{
    float l = cam.near_plane.x;
    float r = cam.near_plane.y;
    float t = cam.near_plane.w;
    float b = cam.near_plane.z;

    parser::Vec3f m = sum_vector(cam.position , multiVector(cam.near_distance,cam.gaze));
    parser::Vec3f v = cam.up;
    parser::Vec3f u = crossVector(v,multiVector(-1,cam.gaze));    
    parser::Vec3f q = sum_vector(sum_vector(m,multiVector(l,u)),multiVector(t,v));

    float su = ((r-l)*( x+0.5))/cam.image_width;
    float sv = ((t-b)*( y+0.5))/cam.image_height;

    parser::Vec3f temp1 = multiVector(su,u);
    parser::Vec3f temp2 = multiVector(sv*(-1),v);
    
    parser::Vec3f s = sum_vector(sum_vector(q,temp1),temp2);
    
    return s;
}

float spheresIntersect(const parser::Vec3f &ray,
                       const parser::Vec3f &pos,
                       const parser::Vec3f &sphereCenter,
                       const float &rad)
{
    parser::Vec3f o_c = subtract_vector(pos,sphereCenter);
    float o_c2 = scalarVector(o_c,o_c);
    float do_c = scalarVector(ray,o_c);
    float do_c2 = do_c*do_c;
    float d2 = scalarVector(ray,ray);
    parser::Vec3f d_neg = multiVector(-1,ray);

    if(do_c2-(d2*(o_c2-rad*rad))>0)
    {
        float sqrt2 = sqrtf(do_c2-(d2*(o_c2-rad*rad)));
        float result1 = scalarVector(d_neg,o_c)-sqrt2;
        float result2 = scalarVector(d_neg,o_c)+sqrt2;
        result1/=d2;
        result2/=d2;

        if(result1>0)
            return result1;
        if(result2>0)
            return result2;
    }
    return std::numeric_limits<float>::max();
}

float determinant(const float &a,
                  const float &b,
                  const float &c,
                  const float &d,
                  const float &e,
                  const float &f,
                  const float &g,
                  const float &h,
                  const float &i)
{
    float result = a*(e*i-h*f)+b*(g*f-d*i)+c*(d*h-e*g);
    return result;
}

float trianglesIntersect(const parser::Vec3f &ray,
                          const parser::Vec3f &position,
                          const std::vector<parser::Vec3f> &vertex_data,
                          const parser::Triangle &triangle)
{
    parser::Vec3f a = vertex_data[triangle.indices.v0_id-1];
    parser::Vec3f b = vertex_data[triangle.indices.v1_id-1];
    parser::Vec3f c = vertex_data[triangle.indices.v2_id-1];

    float ax_ox = a.x-position.x;
    float ay_oy = a.y-position.y;
    float az_oz = a.z-position.z;

    float ax_bx = a.x-b.x;
    float ay_by = a.y-b.y;
    float az_bz = a.z-b.z;

    float ax_cx = a.x-c.x;
    float ay_cy = a.y-c.y;
    float az_cz = a.z-c.z;

    float aDet = determinant(ax_bx,ay_by,az_bz,ax_cx,ay_cy,az_cz,ray.x,ray.y,ray.z);
    float beta = determinant(ax_ox,ay_oy,az_oz,ax_cx,ay_cy,az_cz,ray.x,ray.y,ray.z)/aDet;
    if(beta<0)
    {
        return std::numeric_limits<float>::max();
    }
    float gamma = determinant(ax_bx,ay_by,az_bz,ax_ox,ay_oy,az_oz,ray.x,ray.y,ray.z)/aDet;
    if(gamma<0)
    {
        return std::numeric_limits<float>::max();
    }    
    float t = determinant(ax_bx,ay_by,az_bz,ax_cx,ay_cy,az_cz,ax_ox,ay_oy,az_oz)/aDet;

    if(beta+gamma <= 1 && beta>=0  && gamma >= 0)
    {
        return t;
    }


    return std::numeric_limits<float>::max();


}

float meshHelper(const parser::Vec3f &ray,
                          const parser::Vec3f &position,
                          const std::vector<parser::Vec3f> &vertex_data,
                          const parser::Face &face)
{
    parser::Vec3f a = vertex_data[face.v0_id-1];
    parser::Vec3f b = vertex_data[face.v1_id-1];
    parser::Vec3f c = vertex_data[face.v2_id-1];

    float ax_ox = a.x-position.x;
    float ay_oy = a.y-position.y;
    float az_oz = a.z-position.z;

    float ax_bx = a.x-b.x;
    float ay_by = a.y-b.y;
    float az_bz = a.z-b.z;

    float ax_cx = a.x-c.x;
    float ay_cy = a.y-c.y;
    float az_cz = a.z-c.z;

    float aDet = determinant(ax_bx,ay_by,az_bz,ax_cx,ay_cy,az_cz,ray.x,ray.y,ray.z);
    float beta = determinant(ax_ox,ay_oy,az_oz,ax_cx,ay_cy,az_cz,ray.x,ray.y,ray.z)/aDet;
    if(beta<0)
    {
        return std::numeric_limits<float>::max();
    }
    float gamma = determinant(ax_bx,ay_by,az_bz,ax_ox,ay_oy,az_oz,ray.x,ray.y,ray.z)/aDet;
    if(gamma<0)
    {
        return std::numeric_limits<float>::max();
    }    
    float t = determinant(ax_bx,ay_by,az_bz,ax_cx,ay_cy,az_cz,ax_ox,ay_oy,az_oz)/aDet;

    if(beta+gamma <= 1 && beta>=0  && gamma >= 0)
    {
        return t;
    }

    return std::numeric_limits<float>::max();
}

std::vector<float> meshIntersect(const parser::Vec3f &ray,
                          const parser::Vec3f &position,
                          const std::vector<parser::Vec3f> &vertex_data,
                          const parser::Mesh &mesh)
{
	std::vector<float> v;
    std::vector<parser::Face> faceVector = mesh.faces; 
    float tmin = std::numeric_limits<float>::max();        
    float index;
    for (unsigned int i = 0; i < faceVector.size(); ++i)
    {
        float t = meshHelper(ray,position,vertex_data,faceVector[i]);
        if(t < tmin && t>0)
        {
        	index = i;
            tmin = t ;
        } 
    }
    v.push_back(tmin);
    v.push_back(index);
    return v;

}

std::vector<int> mirrorRecur(const parser::Scene &scene,
                          const parser::Vec3f &cam_pos,
                          const parser::Vec3f &ray,
                          const int depth)
{
    int return_r = 0;
    int return_g = 0;
    int return_b = 0;
    if(depth > scene.max_recursion_depth) {
        std::vector<int> rgb_vector = {scene.background_color.x,scene.background_color.y,scene.background_color.z};
        return rgb_vector;
    }

    float tmin = std::numeric_limits<float>::max();
    float t2min = std::numeric_limits<float>::max();       
    const parser::Sphere * intersectSphere = NULL;
    const parser::Triangle * intersectTriangle = NULL;
    const parser::Mesh * intersectMesh = NULL;

    // SPHERE INTERSECTION
    for (unsigned int k = 0; k<scene.spheres.size(); k++)
    {
        float t = spheresIntersect(ray,cam_pos,scene.vertex_data[scene.spheres[k].center_vertex_id-1],scene.spheres[k].radius);
        if(t<tmin && t > 0)
        {
            tmin = t;
            intersectSphere = &scene.spheres[k];
        }

    }
    // Triangle INTERSECTION
    for (unsigned int k = 0; k<scene.triangles.size(); k++)
    {
        float t = trianglesIntersect(ray,cam_pos,scene.vertex_data,scene.triangles[k]);
        if(t<tmin && t > 0)
        {
            tmin = t ;
            intersectTriangle = &scene.triangles[k];
            intersectSphere = NULL;
        }
    }                
    // Mesh INTERSECTION
    int meshIndex = 0;
    for (unsigned int k = 0; k<scene.meshes.size(); k++)
    {
        std::vector<float> v = meshIntersect(ray,cam_pos,scene.vertex_data,scene.meshes[k]);
        float t = v[0];
        if(t<tmin && t > 0)
        {
            meshIndex = v[1]; 
            tmin = t ;
            intersectMesh = &scene.meshes[k];
            intersectSphere = NULL;
            intersectTriangle = NULL;
        }
    }                                

    float diffuser,diffuseb,diffuseg,specularr,specularb,specularg,mirrorr,mirrorg,mirrorb;
    diffuser = 0; diffuseb = 0; diffuseg = 0; specularr = 0; specularb = 0; specularg = 0; mirrorr = 0 ; mirrorb = 0 ; mirrorg = 0;
    if(intersectSphere)
    {

        float ambientr = scene.materials[intersectSphere->material_id-1].ambient.x* scene.ambient_light.x;
        float ambientg = scene.materials[intersectSphere->material_id-1].ambient.y* scene.ambient_light.y;
        float ambientb = scene.materials[intersectSphere->material_id-1].ambient.z* scene.ambient_light.z;

        parser::Vec3f objectposition = sum_vector(multiVector(tmin,ray),cam_pos);

        parser::Vec3f sphereNormal = divideVector(intersectSphere->radius,subtract_vector(objectposition,scene.vertex_data[intersectSphere->center_vertex_id-1]));
        parser::Vec3f wo =  divideVector(normVector(subtract_vector(cam_pos, objectposition)),subtract_vector(cam_pos, objectposition));

        for (unsigned int indLight = 0; indLight < scene.point_lights.size(); ++indLight)
        {   
            int state = 1;

            parser::Vec3f temp1 = subtract_vector(scene.point_lights[indLight].position,objectposition);
            parser::Vec3f wi = divideVector(normVector(temp1),temp1);
            
            parser::Vec3f wo_wi = sum_vector(wi,wo);

            float normWoWi = normVector(wo_wi);
            parser::Vec3f h =  divideVector(normWoWi,wo_wi);
            
            float norm = normVector(subtract_vector(scene.point_lights[indLight].position,objectposition));
            float norm2 = (norm*norm);
            
            float tmax = scalarVector(subtract_vector(scene.point_lights[indLight].position,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi))),wi)/scalarVector(wi,wi); 
            
            for (unsigned int spheresCount = 0; spheresCount < scene.spheres.size(); ++spheresCount)
            {

                float t = spheresIntersect(wi,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi)),scene.vertex_data[scene.spheres[spheresCount].center_vertex_id-1],scene.spheres[spheresCount].radius);
                if(t<t2min && t>0 && t<tmax)
                {
                    state=0;
                    break;
                }                            
            }
            
            if(state)
            for (unsigned int triangleCount = 0; triangleCount < scene.triangles.size(); ++triangleCount)
            {
                float t = trianglesIntersect(wi,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi)),scene.vertex_data,scene.triangles[triangleCount]);
                if(t<t2min && t>0 && t<tmax)
                {
                    state=0;
                    break;
                }                               
            }
            
            if(state)
            for (unsigned int meshCount = 0; meshCount < scene.meshes.size(); ++meshCount)
            {                    
                std::vector<float> v = meshIntersect(wi,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi)),scene.vertex_data,scene.meshes[meshCount]);
                float t = v[0];
                if(t<t2min && t>0 && t<tmax)
                {
                    state=0;
                    break;
                }                               
            }

            if(state==1)
            {
                float diffScalar = std::max(0.0f,scalarVector(wi,sphereNormal));
                diffuser += (scene.materials[intersectSphere->material_id-1].diffuse.x * diffScalar) * scene.point_lights[indLight].intensity.x/ norm2; 
                diffuseg += (scene.materials[intersectSphere->material_id-1].diffuse.y * diffScalar) * scene.point_lights[indLight].intensity.y/ norm2;
                diffuseb += (scene.materials[intersectSphere->material_id-1].diffuse.z * diffScalar) * scene.point_lights[indLight].intensity.z/ norm2;                        

                float aa = std::max(0.0f,scalarVector(sphereNormal,h));
                float bb = scene.materials[intersectSphere->material_id-1].phong_exponent;
                float specScalar = powf(aa,bb);

                float kk = scene.materials[intersectSphere->material_id-1].specular.x* specScalar;
                float mm = scene.point_lights[indLight].intensity.x;

                specularr +=  kk*mm/norm2;
                specularg += scene.materials[intersectSphere->material_id-1].specular.y* specScalar *scene.point_lights[indLight].intensity.y/norm2;
                specularb += scene.materials[intersectSphere->material_id-1].specular.z* specScalar *scene.point_lights[indLight].intensity.z/norm2;
            }


        }
        if(scene.materials[intersectSphere->material_id-1].mirror.x || scene.materials[intersectSphere->material_id-1].mirror.y || scene.materials[intersectSphere->material_id-1].mirror.z)
        {
            parser::Vec3f sww = multiVector(2, sphereNormal);
            parser::Vec3f wr = sum_vector(multiVector(-1, wo),multiVector(scalarVector(sphereNormal,wo),sww));
            wr = divideVector(normVector(wr), wr);

            std::vector<int> v = mirrorRecur(scene,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,sphereNormal)),wr,depth+1);

            mirrorr += (scene.materials[intersectSphere->material_id-1].mirror.x) *v[0];
            mirrorg += (scene.materials[intersectSphere->material_id-1].mirror.y) *v[1];
            mirrorb += (scene.materials[intersectSphere->material_id-1].mirror.z) *v[2];
        }
        if(ambientr + diffuser + specularr + mirrorr < 255){
            return_r = std::round(ambientr + diffuser+ specularr + mirrorr);
        }
        else{ 
            return_r = 255;
        }
        if (ambientg + diffuseg + specularg + mirrorg < 255){
            return_g = std::round(ambientg + diffuseg+ specularg + mirrorg);
        }
        else{
            return_g = 255;
        }
        if(ambientb + diffuseb + specularb + mirrorb < 255){
            return_b = std::round(ambientb + diffuseb+ specularb + mirrorb);
        }
        else{
            return_b = 255;
        }
    } 

    else if (intersectTriangle)
    {
        float ambientr = scene.materials[intersectTriangle->material_id-1].ambient.x* scene.ambient_light.x;
        float ambientg = scene.materials[intersectTriangle->material_id-1].ambient.y* scene.ambient_light.y;
        float ambientb = scene.materials[intersectTriangle->material_id-1].ambient.z* scene.ambient_light.z;
        
        parser::Vec3f objectposition = sum_vector(multiVector(tmin,ray),cam_pos);
        parser::Vec3f normal = intersectTriangle->indices.normal;
        parser::Vec3f wo =  divideVector(normVector(subtract_vector(cam_pos, objectposition)),subtract_vector(cam_pos, objectposition));


        for (unsigned int indLight = 0; indLight < scene.point_lights.size(); ++indLight)
        {   
            int state = 1;

            parser::Vec3f wi = divideVector(normVector(subtract_vector(scene.point_lights[indLight].position,objectposition)),subtract_vector(scene.point_lights[indLight].position,objectposition));
            parser::Vec3f wo_wi = sum_vector(wi,wo);
            float normWoWi = normVector(wo_wi);
            parser::Vec3f h =  divideVector(normWoWi,wo_wi);
            float norm = normVector(subtract_vector(scene.point_lights[indLight].position,objectposition));
            float norm2 = (norm*norm);
            float tmax = scalarVector(subtract_vector(scene.point_lights[indLight].position,
                                    sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi))),wi)/scalarVector(wi,wi); 
                                 
            for (unsigned int spheresCount = 0; spheresCount < scene.spheres.size(); ++spheresCount)
            {

                float t = spheresIntersect(wi,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi)),scene.vertex_data[scene.spheres[spheresCount].center_vertex_id-1],scene.spheres[spheresCount].radius);
                if(t<t2min && t>0 && t<tmax)
                {
                    state=0;
                    break;
                }                            
            }
            
            if(state)
            for (unsigned int triangleCount = 0; triangleCount < scene.triangles.size(); ++triangleCount)
            {
                float t = trianglesIntersect(wi,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi)),scene.vertex_data,scene.triangles[triangleCount]);
                if(t<t2min && t>0 && t<tmax)
                {
                    state=0;
                    break;
                }                               
            }
            
            if(state)
            for (unsigned int meshCount = 0; meshCount < scene.meshes.size(); ++meshCount)
            {                    
                std::vector<float> v = meshIntersect(wi,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi)),scene.vertex_data,scene.meshes[meshCount]);
                float t = v[0];
                if(t<t2min && t>0 && t<tmax)
                {
                    state=0;
                    break;
                }                               
            }                        

            if(state==1)
            {
                float diffScalar = std::max(0.0f,scalarVector(wi,normal));

                diffuser += (scene.materials[intersectTriangle->material_id-1].diffuse.x * diffScalar) * scene.point_lights[indLight].intensity.x/ norm2; 
                diffuseg += (scene.materials[intersectTriangle->material_id-1].diffuse.y * diffScalar) * scene.point_lights[indLight].intensity.y/ norm2;
                diffuseb += (scene.materials[intersectTriangle->material_id-1].diffuse.z * diffScalar) * scene.point_lights[indLight].intensity.z/ norm2;

                float specScalar = std::pow(std::max(0.0f,scalarVector(normal,h)),scene.materials[intersectTriangle->material_id-1].phong_exponent);

                specularr += (scene.materials[intersectTriangle->material_id-1].specular.x* specScalar) *scene.point_lights[indLight].intensity.x/norm2;
                specularg += (scene.materials[intersectTriangle->material_id-1].specular.y* specScalar) *scene.point_lights[indLight].intensity.y/norm2;
                specularb += (scene.materials[intersectTriangle->material_id-1].specular.z* specScalar) *scene.point_lights[indLight].intensity.z/norm2;
            }


        }
        if(scene.materials[intersectTriangle->material_id-1].mirror.x || scene.materials[intersectTriangle->material_id-1].mirror.y || scene.materials[intersectTriangle->material_id-1].mirror.z)
        {
            parser::Vec3f wr = sum_vector(multiVector(-1,wo),multiVector(2*scalarVector(normal,wo),normal));
            wr = divideVector(normVector(wr), wr);

            std::vector<int> v = mirrorRecur(scene,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,normal)),wr,depth+1);
            mirrorr += (scene.materials[intersectTriangle->material_id-1].mirror.x) *v[0];
            mirrorg += (scene.materials[intersectTriangle->material_id-1].mirror.y) *v[1];
            mirrorb += (scene.materials[intersectTriangle->material_id-1].mirror.z) *v[2];
        }
        if(ambientr + diffuser + specularr + mirrorr < 255){
            return_r = std::round(ambientr + diffuser+ specularr + mirrorr);
        }
        else{ 
            return_r = 255;
        }
        if (ambientg + diffuseg + specularg + mirrorg < 255){
            return_g = std::round(ambientg + diffuseg+ specularg + mirrorg);
        }
        else{
            return_g = 255;
        }
        if(ambientb + diffuseb + specularb + mirrorb < 255){
            return_b = std::round(ambientb + diffuseb+ specularb + mirrorb);
        }
        else{
            return_b = 255;
        }
    }
    else if (intersectMesh)
    {
        float ambientr = scene.materials[intersectMesh->material_id-1].ambient.x* scene.ambient_light.x;
        float ambientg = scene.materials[intersectMesh->material_id-1].ambient.y* scene.ambient_light.y;
        float ambientb = scene.materials[intersectMesh->material_id-1].ambient.z* scene.ambient_light.z;
        
        parser::Vec3f objectposition = sum_vector(multiVector(tmin,ray),cam_pos);
        parser::Vec3f normal = intersectMesh->faces[meshIndex].normal;
        parser::Vec3f wo =  divideVector(normVector(subtract_vector(cam_pos, objectposition)),subtract_vector(cam_pos, objectposition));


        for (unsigned int indLight = 0; indLight < scene.point_lights.size(); ++indLight)
        {
            int state = 1;

            parser::Vec3f wi = divideVector(normVector(subtract_vector(scene.point_lights[indLight].position,objectposition)),
                                            subtract_vector(scene.point_lights[indLight].position,objectposition));
            parser::Vec3f wo_wi = sum_vector(wi,wo);
            float normWoWi = normVector(wo_wi);
            parser::Vec3f h =  divideVector(normWoWi,wo_wi);                        
            float norm = normVector(subtract_vector(scene.point_lights[indLight].position,objectposition));
            float norm2 = (norm*norm);
            float tmax = scalarVector(subtract_vector(scene.point_lights[indLight].position,
                                            sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi))),wi)/scalarVector(wi,wi); 

            for(unsigned int spheresCount = 0; spheresCount < scene.spheres.size(); ++spheresCount)
            {
                float t = spheresIntersect(wi,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi)),scene.vertex_data[scene.spheres[spheresCount].center_vertex_id-1],scene.spheres[spheresCount].radius);
                if(t<t2min && t>0 && t<tmax)
                {
                    state=0;
                    break;
                }                            
            }
            
            if(state)
            for(unsigned int triangleCount = 0; triangleCount < scene.triangles.size(); ++triangleCount)
            {
                float t = trianglesIntersect(wi,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi)),scene.vertex_data,scene.triangles[triangleCount]);
                if(t<t2min && t>0 && t<tmax)
                {
                    state=0;
                    break;
                }                               
            }
            
            if(state)
            for(unsigned int meshCount = 0; meshCount < scene.meshes.size(); ++meshCount)
            {                    
                std::vector<float> v = meshIntersect(wi,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,wi)),scene.vertex_data,scene.meshes[meshCount]);

                float t = v[0];
                if(t<t2min && t>0 && t<tmax)
                {
                    state=0;
                    break;
                }                               
            }                  

            if(state==1){
                float diffScalar = std::max(0.0f,scalarVector(wi,normal));

                diffuser += (scene.materials[intersectMesh->material_id-1].diffuse.x * diffScalar) * scene.point_lights[indLight].intensity.x/ norm2; 
                diffuseg += (scene.materials[intersectMesh->material_id-1].diffuse.y * diffScalar) * scene.point_lights[indLight].intensity.y/ norm2;
                diffuseb += (scene.materials[intersectMesh->material_id-1].diffuse.z * diffScalar) * scene.point_lights[indLight].intensity.z/ norm2; 
                
                float specScalar = std::pow(std::max(0.0f,scalarVector(normal,h)),scene.materials[intersectMesh->material_id-1].phong_exponent);

                specularr += (scene.materials[intersectMesh->material_id-1].specular.x* specScalar) *scene.point_lights[indLight].intensity.x/norm2;
                specularg += (scene.materials[intersectMesh->material_id-1].specular.y* specScalar) *scene.point_lights[indLight].intensity.y/norm2;
                specularb += (scene.materials[intersectMesh->material_id-1].specular.z* specScalar) *scene.point_lights[indLight].intensity.z/norm2; 
            }
                

        }
        if(scene.materials[intersectMesh->material_id-1].mirror.x || scene.materials[intersectMesh->material_id-1].mirror.y || scene.materials[intersectMesh->material_id-1].mirror.z)
        {
            parser::Vec3f wr = sum_vector(multiVector(-1,wo),multiVector(2*scalarVector(normal,wo),normal));
            wr = divideVector(normVector(wr), wr);

            std::vector<int> v = mirrorRecur(scene,sum_vector(objectposition,multiVector(scene.shadow_ray_epsilon,normal)),wr,depth+1);
            mirrorr += (scene.materials[intersectMesh->material_id-1].mirror.x) *v[0];
            mirrorg += (scene.materials[intersectMesh->material_id-1].mirror.y) *v[1];
            mirrorb += (scene.materials[intersectMesh->material_id-1].mirror.z) *v[2];
        }
        if(ambientr + diffuser + specularr + mirrorr < 255){
            return_r = std::round(ambientr + diffuser+ specularr + mirrorr);
        }
        else{ 
            return_r = 255;
        }
        if (ambientg + diffuseg + specularg + mirrorg < 255){
            return_g = std::round(ambientg + diffuseg+ specularg + mirrorg);
        }
        else{
            return_g = 255;
        }
        if(ambientb + diffuseb + specularb + mirrorb < 255){
            return_b = std::round(ambientb + diffuseb+ specularb + mirrorb);
        }
        else{
            return_b = 255;
        }
    }
    else
    {
        return_r = scene.background_color.x;
        return_g = scene.background_color.y;
        return_b = scene.background_color.z;
    }

    std::vector<int> rgb_vector = {return_r, return_g, return_b};
    return rgb_vector;        
}

