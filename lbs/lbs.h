// fix error in old algorithm
#ifndef LBS_H
#define LBS_H

#include <vector>
#include <fstream>
#include "Eigen/Dense"  
using namespace Eigen;
namespace lbs{
    /*
     * Data Types
     */          
    struct Gradient{
        MatrixXf Dx;
        MatrixXf Dy;
    }; 
    struct Operator{
        MatrixXf f2v;
        MatrixXf v2f;
        MatrixXf vertexLaplacianOnVertex;
        Gradient faceGradientOnVertex;
        MatrixXf fullGradientOnVertex;
        MatrixXf fullDivergenceOnFace;
    };  
    // for efficient computation
    struct LsqCache{ 
        MatrixXf lbsMatrix; 
        int vertexNum = 0;
        int faceNum = 0; 
        std::vector<float> areas;
        std::vector<int> vertexIndexAtRectBoundaryHorizontal_Top;
        std::vector<int> vertexIndexAtRectBoundaryHorizontal_Down; 
        std::vector<int> vertexIndexAtRectBoundaryVertical_Left; 
        std::vector<int> vertexIndexAtRectBoundaryVertical_Right; 
    };   
    static LsqCache solveOnRect(float __vertices[][2], int __vertexNum, int __faces[][3], int __faceNum, float __mu[][2], float __bounds[4], float __solution[][2]){
         LsqCache lsqCache;
         lsqCache.vertexNum = __vertexNum;
         lsqCache.faceNum = __faceNum;
        // getAreasAndCounterClockwiseVertex   
        for (int f = 0; f < __faceNum; f++){ 
            float a[2] = {__vertices[__faces[f][0]][0], __vertices[__faces[f][0]][1]};
            float b[2] = {__vertices[__faces[f][1]][0], __vertices[__faces[f][1]][1]};
            float c[2] = {__vertices[__faces[f][2]][0], __vertices[__faces[f][2]][1]};
            float area = 0.5*(a[0]*(b[1]-c[1])+b[0]*(c[1]-a[1])+c[0]*(a[1]-b[1])); 
            if (area < 0){
                 std::swap(__faces[f][0], __faces[f][1]);
                lsqCache.areas.push_back(-area);
            }else{                
                lsqCache.areas.push_back(area);
            } 
        } 
        // get horizontal & vertical boundary vertex    
        lsqCache.vertexIndexAtRectBoundaryHorizontal_Top.push_back(0);
        lsqCache.vertexIndexAtRectBoundaryHorizontal_Down.push_back(0);
        lsqCache.vertexIndexAtRectBoundaryVertical_Left.push_back(0);
        lsqCache.vertexIndexAtRectBoundaryVertical_Right.push_back(0);
        int __tempMinXIndex = 0;
        int __tempMinYIndex = 0;
        int __tempMaxXIndex = 0;
        int __tempMaxYIndex = 0; 
        for (int i = 1; i< __vertexNum; i++){  
            if(__vertices[i][0] < __vertices[__tempMinXIndex][0]){
                __tempMinXIndex = i;
                lsqCache.vertexIndexAtRectBoundaryVertical_Left.clear();
                lsqCache.vertexIndexAtRectBoundaryVertical_Left.push_back(i);
            }else if (__vertices[i][0] == __vertices[__tempMinXIndex][0]){ 
                lsqCache.vertexIndexAtRectBoundaryVertical_Left.push_back(i); 
            }if (__vertices[i][0] > __vertices[__tempMaxXIndex][0]){ 
                __tempMaxXIndex = i;
                lsqCache.vertexIndexAtRectBoundaryVertical_Right.clear();
                lsqCache.vertexIndexAtRectBoundaryVertical_Right.push_back(i);
            }else if (__vertices[i][0] == __vertices[__tempMaxXIndex][0]){  
                lsqCache.vertexIndexAtRectBoundaryVertical_Right.push_back(i);
            }if (__vertices[i][1] > __vertices[__tempMaxYIndex][1]){ 
                __tempMaxYIndex = i;
                lsqCache.vertexIndexAtRectBoundaryHorizontal_Down.clear();
                lsqCache.vertexIndexAtRectBoundaryHorizontal_Down.push_back(i);
            }else if (__vertices[i][1] == __vertices[__tempMaxYIndex][1]){ 
                lsqCache.vertexIndexAtRectBoundaryHorizontal_Down.push_back(i);
            }if(__vertices[i][1] < __vertices[__tempMinYIndex][1]){
                __tempMinYIndex = i;
                lsqCache.vertexIndexAtRectBoundaryHorizontal_Top.clear();
                lsqCache.vertexIndexAtRectBoundaryHorizontal_Top.push_back(i);
            }else if (__vertices[i][1] == __vertices[__tempMinYIndex][1]){ 
                lsqCache.vertexIndexAtRectBoundaryHorizontal_Top.push_back(i); 
            }  
        }       
        float __edges[lsqCache.faceNum][3][2];  
        for(int f=0; f<lsqCache.faceNum;f++){
            __edges[f][0][0] = __vertices[__faces[f][2]][0] - __vertices[__faces[f][1]][0];
            __edges[f][0][1] = __vertices[__faces[f][2]][1] - __vertices[__faces[f][1]][1]; 
            __edges[f][1][0] = __vertices[__faces[f][0]][0] - __vertices[__faces[f][2]][0];
            __edges[f][1][1] = __vertices[__faces[f][0]][1] - __vertices[__faces[f][2]][1]; 
            __edges[f][2][0] = __vertices[__faces[f][1]][0] - __vertices[__faces[f][0]][0];
            __edges[f][2][1] = __vertices[__faces[f][1]][1] - __vertices[__faces[f][0]][1]; 
        }  
        MatrixXf __lbsMatrix(lsqCache.vertexNum, lsqCache.vertexNum);
        __lbsMatrix.setZero();
        for(int f=0; f<lsqCache.faceNum; f++){
            float p = __mu[f][0];
            float r = __mu[f][1];
            float base = 1-p*p-r*r;
            float a_1 = (-(p-1)*(p-1)-r*r)/base;
            float a_2 = 2*r/base; 
            float a_4 = (-(1+p)*(1+p)-r*r)/base;
            for (int i=0; i<3; i++){ 
            for (int j=0; j<3; j++){
                __lbsMatrix(__faces[f][i], __faces[f][j]) += 
                        (a_1*__edges[f][i][1]*__edges[f][j][1] - 
                        a_2*(__edges[f][i][0]*__edges[f][j][1] + __edges[f][i][1]*__edges[f][j][0]) +
                        a_4*__edges[f][i][0]*__edges[f][j][0])/lsqCache.areas[f]; 
            }
            } 
        }   
        lsqCache.lbsMatrix = __lbsMatrix;
        // handle boundary conditions and organize least square vector
        VectorXf __lbsVector_U(__vertexNum);
        VectorXf __lbsVector_V(__vertexNum); 
        __lbsVector_U.setZero();
        __lbsVector_V.setZero(); 
        
        MatrixXf __lbsMatrix_U = lsqCache.lbsMatrix;
        MatrixXf __lbsMatrix_V = lsqCache.lbsMatrix;   
         // horizontal boundary
         float b = __bounds[2];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Top){  
         __lbsVector_V += b* __lbsMatrix_V.col(j);
         for (int i=0; i<__lbsMatrix_V.rows(); i++){  
         __lbsMatrix_V(i,j) = 0;   
         }  
         }b = __bounds[3];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Down){  
         __lbsVector_V += b* __lbsMatrix_V.col(j);
         for (int i=0; i<__lbsMatrix_V.rows(); i++){  
         __lbsMatrix_V(i,j) = 0;   
         } 
         } 
         // vertical boundary
         b = __bounds[0];  
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Left){   
         __lbsVector_U += b* __lbsMatrix_U.col(j); 
         for (int i=0; i<__lbsMatrix_U.rows(); i++){  
         __lbsMatrix_U(i,j) = 0;   
         } 
         }b = __bounds[1];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Right){   
         __lbsVector_U += b* __lbsMatrix_U.col(j);
         for (int i=0; i<__lbsMatrix_U.rows(); i++){  
         __lbsMatrix_U(i,j) = 0;   
         }  
         } 
         // horizontal boundary
         b = __bounds[2];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Top){ 
         for (int i=0; i<__lbsMatrix_V.cols(); i++){  
         __lbsMatrix_V(j,i) = 0;    
         }__lbsMatrix_V(j,j) = 1; 
         __lbsVector_V(j) = -b;
         }b = __bounds[3];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Down){    
         for (int i=0; i<__lbsMatrix_V.cols(); i++){  
         __lbsMatrix_V(j,i) = 0;    
         }__lbsMatrix_V(j,j) = 1;  
         __lbsVector_V(j) = -b;
         } 
         // vertical boundary
         b = __bounds[0];  
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Left){      
         for (int i=0; i<__lbsMatrix_U.cols(); i++){  
         __lbsMatrix_U(j,i) = 0;    
         }__lbsMatrix_U(j,j) = 1;  
         __lbsVector_U(j) = -b;
         }b = __bounds[1];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Right){      
         for (int i=0; i<__lbsMatrix_U.cols(); i++){  
         __lbsMatrix_U(j,i) = 0;    
         }__lbsMatrix_U(j,j) = 1;  
         __lbsVector_U(j) = -b;
         } 
        VectorXf sol_U = __lbsMatrix_U.colPivHouseholderQr().solve(-__lbsVector_U);
        VectorXf sol_V = __lbsMatrix_V.colPivHouseholderQr().solve(-__lbsVector_V);  
        // fill back solution
        for (int i = 0;i< lsqCache.vertexNum;i++){
            __solution[i][0] = sol_U(i);
            __solution[i][1] = sol_V(i);
        } 
         return lsqCache;
    }      
    static LsqCache solveOnRect(LsqCache lsqCache, float __vertices[][2], int __faces[][3], float __mu[][2], float __bounds[4], float __solution[][2]){                
        float __edges[lsqCache.faceNum][3][2];  
        for(int f=0; f<lsqCache.faceNum;f++){
            __edges[f][0][0] = __vertices[__faces[f][2]][0] - __vertices[__faces[f][1]][0];
            __edges[f][0][1] = __vertices[__faces[f][2]][1] - __vertices[__faces[f][1]][1]; 
            __edges[f][1][0] = __vertices[__faces[f][0]][0] - __vertices[__faces[f][2]][0];
            __edges[f][1][1] = __vertices[__faces[f][0]][1] - __vertices[__faces[f][2]][1]; 
            __edges[f][2][0] = __vertices[__faces[f][1]][0] - __vertices[__faces[f][0]][0];
            __edges[f][2][1] = __vertices[__faces[f][1]][1] - __vertices[__faces[f][0]][1]; 
        }   
        lsqCache.lbsMatrix.setZero();
        for(int f=0; f<lsqCache.faceNum; f++){
            float p = __mu[f][0];
            float r = __mu[f][1];
            float base = 1-p*p-r*r;
            float a_1 = (-(p-1)*(p-1)-r*r)/base;
            float a_2 = 2*r/base; 
            float a_4 = (-(1+p)*(1+p)-r*r)/base;
            for (int i=0; i<3; i++){ 
            for (int j=0; j<3; j++){
                lsqCache.lbsMatrix(__faces[f][i], __faces[f][j]) += 
                        (a_1*__edges[f][i][1]*__edges[f][j][1] - 
                        a_2*(__edges[f][i][0]*__edges[f][j][1] + __edges[f][i][1]*__edges[f][j][0]) +
                        a_4*__edges[f][i][0]*__edges[f][j][0])/lsqCache.areas[f]; 
            }
            } 
        }    
        // handle boundary conditions and organize least square vector
        VectorXf __lbsVector_U(lsqCache.vertexNum);
        VectorXf __lbsVector_V(lsqCache.vertexNum); 
        __lbsVector_U.setZero();
        __lbsVector_V.setZero(); 
        
        MatrixXf __lbsMatrix_U = lsqCache.lbsMatrix;
        MatrixXf __lbsMatrix_V = lsqCache.lbsMatrix;   
         // horizontal boundary
         float b = __bounds[2];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Top){  
         __lbsVector_V += b* __lbsMatrix_V.col(j);
         for (int i=0; i<__lbsMatrix_V.rows(); i++){  
         __lbsMatrix_V(i,j) = 0;   
         }  
         }b = __bounds[3];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Down){  
         __lbsVector_V += b* __lbsMatrix_V.col(j);
         for (int i=0; i<__lbsMatrix_V.rows(); i++){  
         __lbsMatrix_V(i,j) = 0;   
         } 
         } 
         // vertical boundary
         b = __bounds[0];  
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Left){   
         __lbsVector_U += b* __lbsMatrix_U.col(j); 
         for (int i=0; i<__lbsMatrix_U.rows(); i++){  
         __lbsMatrix_U(i,j) = 0;   
         } 
         }b = __bounds[1];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Right){   
         __lbsVector_U += b* __lbsMatrix_U.col(j);
         for (int i=0; i<__lbsMatrix_U.rows(); i++){  
         __lbsMatrix_U(i,j) = 0;   
         }  
         } 
         // horizontal boundary
         b = __bounds[2];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Top){ 
         for (int i=0; i<__lbsMatrix_V.cols(); i++){  
         __lbsMatrix_V(j,i) = 0;    
         }__lbsMatrix_V(j,j) = 1; 
         __lbsVector_V(j) = -b;
         }b = __bounds[3];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Down){    
         for (int i=0; i<__lbsMatrix_V.cols(); i++){  
         __lbsMatrix_V(j,i) = 0;    
         }__lbsMatrix_V(j,j) = 1;  
         __lbsVector_V(j) = -b;
         } 
         // vertical boundary
         b = __bounds[0];  
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Left){      
         for (int i=0; i<__lbsMatrix_U.cols(); i++){  
         __lbsMatrix_U(j,i) = 0;    
         }__lbsMatrix_U(j,j) = 1;  
         __lbsVector_U(j) = -b;
         }b = __bounds[1];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Right){      
         for (int i=0; i<__lbsMatrix_U.cols(); i++){  
         __lbsMatrix_U(j,i) = 0;    
         }__lbsMatrix_U(j,j) = 1;  
         __lbsVector_U(j) = -b;
         } 
        VectorXf sol_U = __lbsMatrix_U.colPivHouseholderQr().solve(-__lbsVector_U);
        VectorXf sol_V = __lbsMatrix_V.colPivHouseholderQr().solve(-__lbsVector_V);  
        // fill back solution
        for (int i = 0;i< lsqCache.vertexNum;i++){
            __solution[i][0] = sol_U(i);
            __solution[i][1] = sol_V(i);
        } 
         return lsqCache;
    }

    static LsqCache solveOnRect(LsqCache lsqCache, float __vertices[][2], int __faces[][3], float __bounds[4], float __solution[][2]){       
        
        float __edges[lsqCache.faceNum][3][2];  
        for(int f=0; f<lsqCache.faceNum;f++){
            __edges[f][0][0] = __vertices[__faces[f][2]][0] - __vertices[__faces[f][1]][0];
            __edges[f][0][1] = __vertices[__faces[f][2]][1] - __vertices[__faces[f][1]][1]; 
            __edges[f][1][0] = __vertices[__faces[f][0]][0] - __vertices[__faces[f][2]][0];
            __edges[f][1][1] = __vertices[__faces[f][0]][1] - __vertices[__faces[f][2]][1]; 
            __edges[f][2][0] = __vertices[__faces[f][1]][0] - __vertices[__faces[f][0]][0];
            __edges[f][2][1] = __vertices[__faces[f][1]][1] - __vertices[__faces[f][0]][1]; 
        } 
           
        // handle boundary conditions and organize least square vector
        VectorXf __lbsVector_U(lsqCache.vertexNum);
        VectorXf __lbsVector_V(lsqCache.vertexNum); 
        __lbsVector_U.setZero();
        __lbsVector_V.setZero(); 
        
        MatrixXf __lbsMatrix_U = lsqCache.lbsMatrix;
        MatrixXf __lbsMatrix_V = lsqCache.lbsMatrix;   
         // horizontal boundary
         float b = __bounds[2];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Top){  
         __lbsVector_V += b* __lbsMatrix_V.col(j);
         for (int i=0; i<__lbsMatrix_V.rows(); i++){  
         __lbsMatrix_V(i,j) = 0;   
         }  
         }b = __bounds[3];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Down){  
         __lbsVector_V += b* __lbsMatrix_V.col(j);
         for (int i=0; i<__lbsMatrix_V.rows(); i++){  
         __lbsMatrix_V(i,j) = 0;   
         } 
         } 
         // vertical boundary
         b = __bounds[0];  
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Left){   
         __lbsVector_U += b* __lbsMatrix_U.col(j); 
         for (int i=0; i<__lbsMatrix_U.rows(); i++){  
         __lbsMatrix_U(i,j) = 0;   
         } 
         }b = __bounds[1];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Right){   
         __lbsVector_U += b* __lbsMatrix_U.col(j);
         for (int i=0; i<__lbsMatrix_U.rows(); i++){  
         __lbsMatrix_U(i,j) = 0;   
         }  
         } 
         // horizontal boundary
         b = __bounds[2];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Top){ 
         for (int i=0; i<__lbsMatrix_V.cols(); i++){  
         __lbsMatrix_V(j,i) = 0;    
         }__lbsMatrix_V(j,j) = 1; 
         __lbsVector_V(j) = -b;
         }b = __bounds[3];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Down){    
         for (int i=0; i<__lbsMatrix_V.cols(); i++){  
         __lbsMatrix_V(j,i) = 0;    
         }__lbsMatrix_V(j,j) = 1;  
         __lbsVector_V(j) = -b;
         } 
         // vertical boundary
         b = __bounds[0];  
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Left){      
         for (int i=0; i<__lbsMatrix_U.cols(); i++){  
         __lbsMatrix_U(j,i) = 0;    
         }__lbsMatrix_U(j,j) = 1;  
         __lbsVector_U(j) = -b;
         }b = __bounds[1];
         for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Right){      
         for (int i=0; i<__lbsMatrix_U.cols(); i++){  
         __lbsMatrix_U(j,i) = 0;    
         }__lbsMatrix_U(j,j) = 1;  
         __lbsVector_U(j) = -b;
         } 
        VectorXf sol_U = __lbsMatrix_U.colPivHouseholderQr().solve(-__lbsVector_U);
        VectorXf sol_V = __lbsMatrix_V.colPivHouseholderQr().solve(-__lbsVector_V);  
        // fill back solution
        for (int i = 0;i< lsqCache.vertexNum;i++){
            __solution[i][0] = sol_U(i);
            __solution[i][1] = sol_V(i);
        } 
         return lsqCache;
    }    
    /*
     * Static variables 
     */     
    static Operator operators; 
    static int _vertexNum;
    static int _faceNum;
    static std::vector<std::array<float, 2>> _vertices;
    static std::vector<std::array<int, 3>> _faces;
    static std::vector<float> _areas; 
    static std::vector<std::vector<int>> _facesNearCenterVertex; 
    static std::vector<int> _vertexIndexAtRectBoundaryHorizontal_Top;
    static std::vector<int> _vertexIndexAtRectBoundaryHorizontal_Down; 
    static std::vector<int> _vertexIndexAtRectBoundaryVertical_Left; 
    static std::vector<int> _vertexIndexAtRectBoundaryVertical_Right;  
    /*
     * Methods
     */   
    static void init(float __vertices[][2], int __vertexNum, int __faces[][3], int __faceNum){ 
        // clear variables
        _vertices.clear();
        _faces.clear();
        _facesNearCenterVertex.clear();
        _vertexIndexAtRectBoundaryHorizontal_Top.clear();
        _vertexIndexAtRectBoundaryHorizontal_Down.clear();
        _vertexIndexAtRectBoundaryVertical_Left.clear();
        _vertexIndexAtRectBoundaryVertical_Right.clear(); 
        // re-assign values
        _vertexNum = __vertexNum;
        _faceNum = __faceNum;    
        for (int i=0;i<_vertexNum;i++){           
            std::array<float,2> temp = {__vertices[i][0],__vertices[i][1]};
            _vertices.push_back(temp); 
            std::vector<int> tempVector;
            _facesNearCenterVertex.push_back(tempVector); 
        }for (int i=0;i<_faceNum;i++){        
            std::array<int,3> temp = {__faces[i][0],__faces[i][1],__faces[i][2]}; 
            _faces.push_back(temp); 
        }  
    }
    static void getAreasAndCounterClockwiseVertex(){
        for (int f = 0; f < _faceNum; f++){ 
            float a[2] = {_vertices[_faces[f][0]][0], _vertices[_faces[f][0]][1]};
            float b[2] = {_vertices[_faces[f][1]][0], _vertices[_faces[f][1]][1]};
            float c[2] = {_vertices[_faces[f][2]][0], _vertices[_faces[f][2]][1]};
            float area = 0.5*(a[0]*(b[1]-c[1])+b[0]*(c[1]-a[1])+c[0]*(a[1]-b[1])); 
            if (area < 0){
                 std::swap(_faces[f][0], _faces[f][1]);
                _areas.push_back(-area);
            }else{   
                _areas.push_back(area);
            } 
        } 
    }
    static void createFaceGradient(){  
        operators.faceGradientOnVertex.Dx = MatrixXf(_faceNum, _vertexNum);
        operators.faceGradientOnVertex.Dy = MatrixXf(_faceNum, _vertexNum);
        operators.faceGradientOnVertex.Dx.setZero();
        operators.faceGradientOnVertex.Dy.setZero(); 
        for (int f = 0;f < _faceNum;f++){ 
            int i = _faces[f][0];
            int j = _faces[f][1];
            int k = _faces[f][2];  
            float a[2] = {_vertices[j][0] - _vertices[i][0], _vertices[j][1] - _vertices[i][1]};
            float b[2] = {_vertices[k][0] - _vertices[j][0], _vertices[k][1] - _vertices[j][1]};
            float c[2] = {_vertices[i][0] - _vertices[k][0], _vertices[i][1] - _vertices[k][1]};  
            operators.faceGradientOnVertex.Dx(f,i) = -b[1]*0.5/_areas[f];
            operators.faceGradientOnVertex.Dx(f,j) = -c[1]*0.5/_areas[f];
            operators.faceGradientOnVertex.Dx(f,k) = -a[1]*0.5/_areas[f];
            operators.faceGradientOnVertex.Dy(f,i) = b[0]*0.5/_areas[f];
            operators.faceGradientOnVertex.Dy(f,j) = c[0]*0.5/_areas[f];
            operators.faceGradientOnVertex.Dy(f,k) = a[0]*0.5/_areas[f];  
        }   
    }  
    static void createFullGradient(){ 
        operators.fullGradientOnVertex = MatrixXf(_faceNum*2, _vertexNum); 
        operators.fullGradientOnVertex << operators.faceGradientOnVertex.Dx, operators.faceGradientOnVertex.Dy;
    }
    static void createFullDivergence(){ 
        operators.fullDivergenceOnFace = MatrixXf(_vertexNum, _faceNum*2);
        operators.fullDivergenceOnFace = operators.fullGradientOnVertex.transpose();
    }
    static float dot(float v1[2], float v2[2]){
        return v1[0]*v2[0] + v1[1]*v2[1];
    }
    static float cotangent(std::array<float, 2> centerVertex, std::array<float, 2> angleVertex, std::array<float, 2> outerVertex){
        float a[2],b[2];
        a[0] = centerVertex[0] - angleVertex[0];
        a[1] = centerVertex[1] - angleVertex[1];
        b[0] = outerVertex[0] - angleVertex[0];
        b[1] = outerVertex[1] - angleVertex[1];
        return dot(a, b)/std::abs(a[0]*b[1] - a[1]*b[0]);
    }
    static void vertexLaplacian(){
        float laplacianArea[_vertexNum];
        std::fill(laplacianArea, laplacianArea + _vertexNum, 0); 
        operators.vertexLaplacianOnVertex = MatrixXf(_vertexNum, _vertexNum);
        for (int f = 0; f<_faceNum;f++){
            int centre = 0;
            int outer = 1;
            int angle = 2; 
            
            for (int t = 0; t<3; t++){     
                laplacianArea[_faces[f][centre]] += _areas[f]/4; 
                for (int tt = 0; tt<2; tt++){ 
                    float temp = cotangent(_vertices[_faces[f][centre]], _vertices[_faces[f][angle]], _vertices[_faces[f][outer]]); 
                    operators.vertexLaplacianOnVertex(_faces[f][centre], _faces[f][centre]) += temp;
                    operators.vertexLaplacianOnVertex(_faces[f][centre], _faces[f][outer]) -= temp;  
                    std::swap(outer, angle);
                } 
                std::swap(outer, angle); 
                centre++; centre%=3;
                outer++; outer%=3;
                angle++; angle%=3; 
            }  
        } 
        // divide back areas
        for (int i = 0; i < _vertexNum; i++){ 
            for (int j = 0; j < _vertexNum; j++){ 
                operators.vertexLaplacianOnVertex(i,j) /= laplacianArea[i];
            }  
        } 
    }
    static void createF2V(){  
        operators.f2v = MatrixXf(_vertexNum, _faceNum);
        for (int f = 0; f < _faceNum; f++){
            _facesNearCenterVertex[_faces[f][0]].push_back(f);
            _facesNearCenterVertex[_faces[f][1]].push_back(f);
            _facesNearCenterVertex[_faces[f][2]].push_back(f);
        }
        for (int v = 0; v < _vertexNum; v++){
            float __faceAreasAroundVertex = 0;
            for (int i = 0; i < _facesNearCenterVertex[v].size(); i++){
                __faceAreasAroundVertex += _areas[_facesNearCenterVertex[v][i]];
            }
            for (int i = 0; i < _facesNearCenterVertex[v].size(); i++){
                operators.f2v(v,_facesNearCenterVertex[v][i]) = _areas[_facesNearCenterVertex[v][i]]/__faceAreasAroundVertex;
            }
        }
    }
    static void createV2F(){
        operators.v2f = MatrixXf(_faceNum, _vertexNum); 
        float __weight = 1.0f/3;
        for (int f = 0; f < _faceNum; f++){
            operators.v2f(f, _faces[f][0]) = __weight;
            operators.v2f(f, _faces[f][1]) = __weight;
            operators.v2f(f, _faces[f][2]) = __weight;
        }
    }   
    
    // get  boundary vertex 
    static Operator createOperators(float __vertices[][2], int __vertexNum, int __faces[][3], int __faceNum){   
        init(__vertices, __vertexNum, __faces, __faceNum);
        getAreasAndCounterClockwiseVertex(); 
        createFaceGradient(); 
        createFullGradient();
        createFullDivergence();
        createF2V();
        createV2F(); 
        vertexLaplacian();   
        return operators;
    };  
    /* edit on 8-4-2016 */
    /*
    static void matrix2Array(MatrixXf A, float* array2d){
        for (int i = 0; i<A.rows(); i++){
            for (int j = 0; j<A.cols(); j++){
                array2d[i*A.cols() + j] = A(i,j);
            }
        }
    }   
    static void vector2Array(VectorXf V, float* array2d){
        for (int i = 0; i<V.rows(); i++){
            array2d[i] = V(i);
        }
    }  
    static void array2Vector(float* array2d, VectorXf V){
        for (int i = 0; i<V.rows(); i++){
            V(i) = array2d[i];
        }
    }  
    static void array2Matrix(float* array2d, MatrixXf& A){
        for (int i = 0; i<A.rows(); i++){
            for (int j = 0; j<A.cols(); j++){
                 A(i,j) = array2d[i*A.cols() + j];
            }
        }
    }  
    
    static void solveLinearSystem(float* _2dArrayMatrixPtr, float _2dArrayVector[], int rows, int cols, float solution[]){
        MatrixXf A(rows, cols);
        array2Matrix(_2dArrayMatrixPtr, A);
        VectorXf b(rows); 
        array2Vector(_2dArrayVector, b);
        vector2Array(A.colPivHouseholderQr().solve(-__lbsVector_U), solution);
    }
     * */
    // bound assignment: [x_min][x_max][y_min][y_max]  
} 
  
#endif /* LBS_H */
