#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "helper.h"
#include <cmath>
#include <limits>
#include <thread>
#include <mutex>

std::mutex m;

typedef unsigned char RGB[3];
int total_y = 0;

void calculate_row(const int &width,const parser::Camera &cam,const parser::Scene &scene, unsigned char* image){

    while(1) {
        int exit = 0;
        int thread_y=0;
        m.lock();
        if(total_y>=0){
            thread_y = total_y;
            total_y--;
        } else {
            exit = 1;
        }
        m.unlock();
        if(thread_y<0 || exit)
            break;
        for (int x = 0; x < width; ++x) {
            parser::Vec3f ray = computeRayEquation(thread_y, x, cam);
            std::vector<int> result_rgb = mirrorRecur(scene, cam.position, ray, 0);
            image[thread_y * width * 3 + x * 3] = result_rgb[0];
            image[thread_y * width * 3 + x * 3 + 1] = result_rgb[1];
            image[thread_y * width * 3 + x * 3 + 2] = result_rgb[2];
        }
    }
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    for (unsigned int i = 0 ;i < scene.cameras.size() ; i++ )
    {   
        parser::Camera cam = scene.cameras[i];
        int width = cam.image_width;
        int height = cam.image_height;
        unsigned char* image = new unsigned char [width * height * 3];

        total_y = height-1;

        int total_thread = 8;

        std::thread *th = new std::thread[total_thread];
        for (int y = 0; y < total_thread; ++y)
        {
            th[y] = std::thread(&calculate_row, width, cam, scene, std::ref(image));
//            calculate_row(y,width,cam,scene,image);
        }

        for (int k = 0; k < total_thread; ++k) {
            th[k].join();
        }

        const char* img_name = cam.image_name.c_str();

        write_ppm(img_name, image, width, height);
        delete [] image;
        delete [] th;
    }
    
}