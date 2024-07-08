#include "setGPU.h"
#include <cuda.h>
#include <iostream>

#include <vector>
#include <boost/atomic.hpp>
static boost::atomic<int> activeCudaDeviceIndex;
static std::vector<int> availableGPUs;
bool initGPUSelection()
{
    int devicesCount;
    cudaGetDeviceCount(&devicesCount);
    for(int deviceIndex = 0; deviceIndex < devicesCount; ++deviceIndex)
    {
        cudaDeviceProp deviceProperties;
        cudaGetDeviceProperties(&deviceProperties, deviceIndex);
        std::cout << "Found "<< deviceProperties.name << " PCI id " << deviceProperties.pciBusID << std::endl;
        if (deviceProperties.major >= 3
            && deviceProperties.minor >= 0)
        {
            std::cout << "Adding "<< deviceProperties.name << " PCI id " << deviceProperties.pciBusID << " to available GPUs" << std::endl;
            availableGPUs.push_back(deviceIndex);
            if(deviceProperties.major < 6){
                std::cout << "Newer GPU architecture recommended (6+)" << std::endl;
            }
        }
    }

    if(availableGPUs.size() == 0){
    std::cerr << "New GPU not found! Exiting." << std::endl;
    exit(1);
    return false;
    }
    else{
        activeCudaDeviceIndex = 0;
        return true;
    }
}

void setGPU(int GPUid){
    cudaSetDevice(GPUid);
}

int getGPU(){//int particleNumber, int pCaNum){
    int deviceIndex = availableGPUs[activeCudaDeviceIndex++ % availableGPUs.size()];
    //cudaDeviceProp deviceProperties;
    //cudaGetDeviceProperties(&deviceProperties, deviceIndex);
    //std::cout << "Selecting " << deviceProperties.name << " PCI id " << deviceProperties.pciBusID << " for particle number " << particleNumber << " and cc number " << pCaNum << std::endl;
    return deviceIndex;
}

void printCurrentGPU(){
    int device;
    cudaGetDevice(&device);

    cudaDeviceProp deviceProperties;
    cudaGetDeviceProperties(&deviceProperties, device);
    std::cout << "Current GPU is "<< deviceProperties.name << " PCI id " << deviceProperties.pciBusID << std::endl;
}
