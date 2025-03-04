#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include <cmath>

using namespace ROOT::Math;

int Vertex(){
    
    PositionVector3D<Cartesian3D<double>> pv(0, 0, 0); 
    PositionVector3D<Cartesian3D<double>> sv(1, 1, 1);  
    DisplacementVector3D<Cartesian3D<double>> p_MuMu(1, 1, 1);  

    double s = ((pv - sv).Dot(p_MuMu)) / p_MuMu.Mag2();


    PositionVector3D<Cartesian3D<double>> pca = sv;
    DisplacementVector3D<Cartesian3D<double>> y = s*p_MuMu;
    pca += y;  

    double minDistance = sqrt((pca - pv).Mag2());

    std::cout << "PV: (" << pv.X() << ", " << pv.Y() << ", " << pv.Z() << ")" << std::endl;
    std::cout << "PCA: (" << pca.X() << ", " << pca.Y() << ", " << pca.Z() << ")" << std::endl;
    std::cout << "Minimalna odległość: " << minDistance << std::endl;

    return 0;
}