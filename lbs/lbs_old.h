#include <vector>
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
        MatrixXf laplacianOnFace;
        Gradient gradientOnVertex;
        MatrixXf fullGradientOnVertex;
        MatrixXf fullDivergenceOnFace;
    };  
    // for efficient computation
    struct LsqCache{
        bool hasCache = false; 
        int vertexNum = 0;
        int faceNum = 0;
        MatrixXf lsqMatrix_U;
        MatrixXf lsqMatrix_V; 
        MatrixXf fullGradientOnVertex;
        MatrixXf fullDivergenceOnFace;
        std::vector<int> vertexIndexAtRectBoundaryHorizontal_Top;
        std::vector<int> vertexIndexAtRectBoundaryHorizontal_Down; 
        std::vector<int> vertexIndexAtRectBoundaryVertical_Left; 
        std::vector<int> vertexIndexAtRectBoundaryVertical_Right; 
    }; 
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
        operators.gradientOnVertex.Dx = MatrixXf(_faceNum, _vertexNum);
        operators.gradientOnVertex.Dy = MatrixXf(_faceNum, _vertexNum);
        for (int f = 0;f < _faceNum;f++){ 
            int i = _faces[f][0];
            int j = _faces[f][1];
            int k = _faces[f][2];  
            float a[2] = {_vertices[j][0] - _vertices[i][0], _vertices[j][1] - _vertices[i][1]};
            float b[2] = {_vertices[k][0] - _vertices[j][0], _vertices[k][1] - _vertices[j][1]};
            float c[2] = {_vertices[i][0] - _vertices[k][0], _vertices[i][1] - _vertices[k][1]};   
            operators.gradientOnVertex.Dx(f,i) = b[1]*0.5f/_areas[f];
            operators.gradientOnVertex.Dx(f,j) = c[1]*0.5f/_areas[f];
            operators.gradientOnVertex.Dx(f,k) = a[1]*0.5f/_areas[f];
            operators.gradientOnVertex.Dy(f,i) = -b[0]*0.5f/_areas[f];
            operators.gradientOnVertex.Dy(f,j) = -c[0]*0.5f/_areas[f];
            operators.gradientOnVertex.Dy(f,k) = -a[0]*0.5f/_areas[f]; 
        }   
    }  
    static void createFullGradient(){ 
        operators.fullGradientOnVertex = MatrixXf(_faceNum*2, _vertexNum);
        operators.fullGradientOnVertex << operators.gradientOnVertex.Dx, operators.gradientOnVertex.Dy;
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
        operators.laplacianOnFace = MatrixXf(_vertexNum, _vertexNum);
        for (int f = 0; f<_faceNum;f++){
            int centre = 0;
            int outer = 1;
            int angle = 2; 
            
            for (int t = 0; t<3; t++){     
                laplacianArea[_faces[f][centre]] += _areas[f]/4; 
                for (int tt = 0; tt<2; tt++){ 
                    float temp = cotangent(_vertices[_faces[f][centre]], _vertices[_faces[f][angle]], _vertices[_faces[f][outer]]); 
                    operators.laplacianOnFace(_faces[f][centre], _faces[f][centre]) += temp;
                    operators.laplacianOnFace(_faces[f][centre], _faces[f][outer]) -= temp;  
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
                operators.laplacianOnFace(i,j) /= laplacianArea[i];
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
    static void getRectBoundaryVertex(){ 
        // get horizontal & vertical boundary vertex   
        _vertexIndexAtRectBoundaryVertical_Left.push_back(0);
        _vertexIndexAtRectBoundaryVertical_Right.push_back(0);
        _vertexIndexAtRectBoundaryHorizontal_Top.push_back(0);
        _vertexIndexAtRectBoundaryHorizontal_Down.push_back(0);
        int __tempMinXIndex = 0;
        int __tempMinYIndex = 0;
        int __tempMaxXIndex = 0;
        int __tempMaxYIndex = 0; 
        for (int i = 1; i< _vertexNum; i++){  
            if(_vertices[i][0] < _vertices[__tempMinXIndex][0]){
                __tempMinXIndex = i;
                _vertexIndexAtRectBoundaryVertical_Left.clear();
                _vertexIndexAtRectBoundaryVertical_Left.push_back(i);
            }else if (_vertices[i][0] == _vertices[__tempMinXIndex][0]){ 
                _vertexIndexAtRectBoundaryVertical_Left.push_back(i); 
            }if (_vertices[i][0] > _vertices[__tempMaxXIndex][0]){ 
                __tempMaxXIndex = i;
                _vertexIndexAtRectBoundaryVertical_Right.clear();
                _vertexIndexAtRectBoundaryVertical_Right.push_back(i);
            }else if (_vertices[i][0] == _vertices[__tempMaxXIndex][0]){  
                _vertexIndexAtRectBoundaryVertical_Right.push_back(i);
            }if (_vertices[i][1] > _vertices[__tempMaxYIndex][1]){ 
                __tempMaxYIndex = i;
                _vertexIndexAtRectBoundaryHorizontal_Down.clear();
                _vertexIndexAtRectBoundaryHorizontal_Down.push_back(i);
            }else if (_vertices[i][1] == _vertices[__tempMaxYIndex][1]){ 
                _vertexIndexAtRectBoundaryHorizontal_Down.push_back(i);
            }if(_vertices[i][1] < _vertices[__tempMinYIndex][1]){
                __tempMinYIndex = i;
                _vertexIndexAtRectBoundaryHorizontal_Top.clear();
                _vertexIndexAtRectBoundaryHorizontal_Top.push_back(i);
            }else if (_vertices[i][1] == _vertices[__tempMinYIndex][1]){ 
                _vertexIndexAtRectBoundaryHorizontal_Top.push_back(i); 
            }  
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
    static void matrix2Array(MatrixXf A, float* array2d){
        for (int i = 0; i<A.rows(); i++){
            for (int j = 0; j<A.cols(); j++){
                array2d[i + j*A.cols()] = A(i,j);
            }
        }
    }
    static void restoreBoundaryValues(LsqCache lsqCache, float __bounds[4], float __solution[][2]){ 
        // horizontal boundary
        for (int &i : lsqCache.vertexIndexAtRectBoundaryHorizontal_Top){  
            __solution[i][1] =  __bounds[2]; 
        }for (int &i : lsqCache.vertexIndexAtRectBoundaryHorizontal_Down){ 
            __solution[i][1] =  __bounds[3]; 
        } 
        // vertical boundary
        for (int &i : lsqCache.vertexIndexAtRectBoundaryVertical_Left){    
            __solution[i][0] =  __bounds[0]; 
        }for (int &i : lsqCache.vertexIndexAtRectBoundaryVertical_Right){    
            __solution[i][0] =  __bounds[1]; 
        }   
    }
    static void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
    {
        unsigned int numRows = matrix.rows();
        unsigned int numCols = matrix.cols()-1;
        if( colToRemove < numCols )
            matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);
        matrix.conservativeResize(numRows,numCols);
    }
    static LsqCache solveOnRect(float __vertices[][2], int __vertexNum, int __faces[][3], int __faceNum, float __mu[][2], float __bounds[4], float __solution[][2]){
        LsqCache lsqCache;  
        init(__vertices, __vertexNum, __faces, __faceNum); 
        lsqCache.faceNum = _faceNum;
        lsqCache.vertexNum = _vertexNum;
        getAreasAndCounterClockwiseVertex(); 
        createFaceGradient();
        createFullGradient();
        createFullDivergence();
        getRectBoundaryVertex();
        lsqCache.fullGradientOnVertex = operators.fullGradientOnVertex;
        lsqCache.fullDivergenceOnFace = operators.fullDivergenceOnFace;
        lsqCache.vertexIndexAtRectBoundaryHorizontal_Top = _vertexIndexAtRectBoundaryHorizontal_Top;
        lsqCache.vertexIndexAtRectBoundaryHorizontal_Down = _vertexIndexAtRectBoundaryHorizontal_Down;
        lsqCache.vertexIndexAtRectBoundaryVertical_Left = _vertexIndexAtRectBoundaryVertical_Left;
        lsqCache.vertexIndexAtRectBoundaryVertical_Right = _vertexIndexAtRectBoundaryVertical_Right;
        MatrixXf __muMatrix_U(lsqCache.faceNum* 2, lsqCache.faceNum* 2); 
        MatrixXf __muMatrix_V(lsqCache.faceNum* 2, lsqCache.faceNum* 2);
        __muMatrix_U.setZero();
        __muMatrix_V.setZero(); 
        // construct mu matrix from mu, assume not exist |mu| >= 1
        for (int face= 0; face< lsqCache.faceNum; face++){
            float p = __mu[face][0];
            float r = __mu[face][1];
            float base = 1-p*p-r*r;
            float a_1 = (-(p-1)*(p-1)-r*r)/base;
            float a_2 = 2*r/base;
            float a_3 = a_2;   
            float a_4 = (-(1+p)*(1+p)-r*r)/base;
            // as they are in negative symmetry, the inverse is negative of the original  
            __muMatrix_U(face, face) = a_1;
            __muMatrix_U(face, lsqCache.faceNum+face) = a_2;
            __muMatrix_U(lsqCache.faceNum+face, face) = a_3;
            __muMatrix_U(lsqCache.faceNum+face, lsqCache.faceNum+face) = a_4;
            __muMatrix_V(face, face) = -a_1; 
            __muMatrix_V(face, lsqCache.faceNum+face) = -a_2;
            __muMatrix_V(lsqCache.faceNum+face, face) = -a_3;
            __muMatrix_V(lsqCache.faceNum+face, lsqCache.faceNum+face) = -a_4; 
        }  
        // compose final least square matrix
        lsqCache.lsqMatrix_U = lsqCache.fullDivergenceOnFace* __muMatrix_U* lsqCache.fullGradientOnVertex;
        lsqCache.lsqMatrix_V = lsqCache.fullDivergenceOnFace* __muMatrix_V* lsqCache.fullGradientOnVertex;
        lsqCache.hasCache = true; 
        MatrixXf __lsqMatrix_U_temp = lsqCache.lsqMatrix_U;
        MatrixXf __lsqMatrix_V_temp = lsqCache.lsqMatrix_V;
        // handle boundary conditions and organize least square vector
        VectorXf __lsqVector_U(lsqCache.vertexNum);
        VectorXf __lsqVector_V(lsqCache.vertexNum); 
        __lsqVector_U.setZero();
        __lsqVector_V.setZero();
        // horizontal boundary
        float b = __bounds[2];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Top){   
            __lsqVector_V += b* __lsqMatrix_V_temp.col(j);
            for (int i=0; i<__lsqMatrix_V_temp.rows(); i++){  
                    __lsqMatrix_V_temp(i,j) = 0;   
            }  
        }b = __bounds[3];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Down){    
            __lsqVector_V += b* __lsqMatrix_V_temp.col(j);
            for (int i=0; i<__lsqMatrix_V_temp.rows(); i++){  
                    __lsqMatrix_V_temp(i,j) = 0;   
            }  
        } 
        // vertical boundary
        b = __bounds[0];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Left){   
            __lsqVector_U += b* __lsqMatrix_U_temp.col(j);
            for (int i=0; i<__lsqMatrix_U_temp.rows(); i++){  
                    __lsqMatrix_U_temp(i,j) = 0;   
            }  
        }b = __bounds[1];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Right){  
            __lsqVector_U += b* __lsqMatrix_U_temp.col(j);
            for (int i=0; i<__lsqMatrix_U_temp.rows(); i++){  
                    __lsqMatrix_U_temp(i,j) = 0;   
            }  
        }   
        VectorXf sol_U = __lsqMatrix_U_temp.jacobiSvd(ComputeThinU | ComputeThinV).solve(-__lsqVector_U);
        VectorXf sol_V = __lsqMatrix_V_temp.jacobiSvd(ComputeThinU | ComputeThinV).solve(-__lsqVector_V); 
        // std::cout << sol_U.format(CleanFmt) << sep << std::endl; 
        // fill back solution
        for (int i = 0;i< lsqCache.vertexNum;i++){
            __solution[i][0] = sol_U(i);
            __solution[i][1] = sol_V(i);
        }
        // fill back boundary conditions   
        restoreBoundaryValues(lsqCache, __bounds, __solution);
        return lsqCache;
    }
    static LsqCache solveOnRect(LsqCache lsqCache, float __mu[][2], float __bounds[4], float __solution[][2]){
        if (lsqCache.hasCache){ 
        MatrixXf __muMatrix_U(lsqCache.faceNum* 2, lsqCache.faceNum* 2); 
        MatrixXf __muMatrix_V(lsqCache.faceNum* 2, lsqCache.faceNum* 2);
        __muMatrix_U.setZero();
        __muMatrix_V.setZero(); 
        // construct mu matrix from mu, assume not exist |mu| >= 1
        for (int face= 0; face< lsqCache.faceNum; face++){
            float p = __mu[face][0];
            float r = __mu[face][1];
            float base = 1-p*p-r*r;
            float a_1 = (-(p-1)*(p-1)-r*r)/base;
            float a_2 = 2*r/base;
            float a_3 = a_2;   
            float a_4 = (-(1+p)*(1+p)-r*r)/base;
            // as they are in negative symmetry, the inverse is negative of the original  
            __muMatrix_U(face, face) = a_1;
            __muMatrix_U(face, lsqCache.faceNum+face) = a_2;
            __muMatrix_U(lsqCache.faceNum+face, face) = a_3;
            __muMatrix_U(lsqCache.faceNum+face, lsqCache.faceNum+face) = a_4;
            __muMatrix_V(face, face) = -a_1; 
            __muMatrix_V(face, lsqCache.faceNum+face) = -a_2;
            __muMatrix_V(lsqCache.faceNum+face, face) = -a_3;
            __muMatrix_V(lsqCache.faceNum+face, lsqCache.faceNum+face) = -a_4; 
        }  
        // compose final least square matrix
        lsqCache.lsqMatrix_U = lsqCache.fullDivergenceOnFace* __muMatrix_U* lsqCache.fullGradientOnVertex;
        lsqCache.lsqMatrix_V = lsqCache.fullDivergenceOnFace* __muMatrix_V* lsqCache.fullGradientOnVertex;
        lsqCache.hasCache = true;  
        MatrixXf __lsqMatrix_U_temp = lsqCache.lsqMatrix_U;
        MatrixXf __lsqMatrix_V_temp = lsqCache.lsqMatrix_V;
        // handle boundary conditions and organize least square vector
        VectorXf __lsqVector_U(lsqCache.vertexNum);
        VectorXf __lsqVector_V(lsqCache.vertexNum); 
        __lsqVector_U.setZero();
        __lsqVector_V.setZero();
        // horizontal boundary
        float b = __bounds[2];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Top){   
            __lsqVector_V += b* __lsqMatrix_V_temp.col(j);
            for (int i=0; i<__lsqMatrix_V_temp.rows(); i++){  
                    __lsqMatrix_V_temp(i,j) = 0;   
            }  
        }b = __bounds[3];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Down){    
            __lsqVector_V += b* __lsqMatrix_V_temp.col(j);
            for (int i=0; i<__lsqMatrix_V_temp.rows(); i++){  
                    __lsqMatrix_V_temp(i,j) = 0;   
            }  
        } 
        // vertical boundary
        b = __bounds[0];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Left){   
            __lsqVector_U += b* __lsqMatrix_U_temp.col(j);
            for (int i=0; i<__lsqMatrix_U_temp.rows(); i++){  
                    __lsqMatrix_U_temp(i,j) = 0;   
            }  
        }b = __bounds[1];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Right){  
            __lsqVector_U += b* __lsqMatrix_U_temp.col(j);
            for (int i=0; i<__lsqMatrix_U_temp.rows(); i++){  
                    __lsqMatrix_U_temp(i,j) = 0;   
            }  
        }  
        VectorXf sol_U = __lsqMatrix_U_temp.jacobiSvd(ComputeThinU | ComputeThinV).solve(-__lsqVector_U);
        VectorXf sol_V = __lsqMatrix_V_temp.jacobiSvd(ComputeThinU | ComputeThinV).solve(-__lsqVector_V);  
            // fill back solution
            for (int i = 0;i< lsqCache.vertexNum;i++){
                __solution[i][0] = sol_U(i);
                __solution[i][1] = sol_V(i);
            }
            // fill back boundary conditions   
            restoreBoundaryValues(lsqCache, __bounds, __solution);
        }return lsqCache; 
    }
    
    static LsqCache solveOnRect(LsqCache lsqCache, float __bounds[4], float __solution[][2]){ 
        if (lsqCache.hasCache){  
            MatrixXf __lsqMatrix_U_temp = lsqCache.lsqMatrix_U;
            MatrixXf __lsqMatrix_V_temp = lsqCache.lsqMatrix_V;
            // handle boundary conditions and organize least square vector
            VectorXf __lsqVector_U(lsqCache.vertexNum);
            VectorXf __lsqVector_V(lsqCache.vertexNum); 
            __lsqVector_U.setZero();
            __lsqVector_V.setZero(); 
        // horizontal boundary
        float b = __bounds[2];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Top){   
            __lsqVector_V += b* __lsqMatrix_V_temp.col(j);
            for (int i=0; i<__lsqMatrix_V_temp.rows(); i++){  
                    __lsqMatrix_V_temp(i,j) = 0;   
            }  
        }b = __bounds[3];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryHorizontal_Down){    
            __lsqVector_V += b* __lsqMatrix_V_temp.col(j);
            for (int i=0; i<__lsqMatrix_V_temp.rows(); i++){  
                    __lsqMatrix_V_temp(i,j) = 0;   
            }  
        } 
        // vertical boundary
        b = __bounds[0];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Left){   
            __lsqVector_U += b* __lsqMatrix_U_temp.col(j);
            for (int i=0; i<__lsqMatrix_U_temp.rows(); i++){  
                    __lsqMatrix_U_temp(i,j) = 0;   
            }  
        }b = __bounds[1];
        for (int &j : lsqCache.vertexIndexAtRectBoundaryVertical_Right){  
            __lsqVector_U += b* __lsqMatrix_U_temp.col(j);
            for (int i=0; i<__lsqMatrix_U_temp.rows(); i++){  
                    __lsqMatrix_U_temp(i,j) = 0;   
            }  
        }  
        VectorXf sol_U = __lsqMatrix_U_temp.jacobiSvd(ComputeThinU | ComputeThinV).solve(-__lsqVector_U);
        VectorXf sol_V = __lsqMatrix_V_temp.jacobiSvd(ComputeThinU | ComputeThinV).solve(-__lsqVector_V);             
            // fill back solution
            for (int i = 0;i< lsqCache.vertexNum;i++){
                __solution[i][0] = sol_U(i);
                __solution[i][1] = sol_V(i);
            }
            // fill back boundary conditions   
            restoreBoundaryValues(lsqCache, __bounds, __solution);
        }return lsqCache;
    }  
    // this function is optional, can be used after getting operators
    // bound assignment: [x_min][x_max][y_min][y_max]
    static void solveOnRect(float __mu[][2], float __bounds[4], float __solution[][2]){  
        MatrixXf __lsqMatrix_U(_vertexNum, _vertexNum); 
        MatrixXf __lsqMatrix_V(_vertexNum, _vertexNum); 
        __lsqMatrix_U.setZero();
        __lsqMatrix_V.setZero();
        MatrixXf __muMatrix_U(_faceNum* 2, _faceNum* 2); 
        MatrixXf __muMatrix_V(_faceNum* 2, _faceNum* 2);
        __muMatrix_U.setZero();
        __muMatrix_V.setZero();
        // construct mu matrix from mu, assume not exist |mu| >= 1
        for (int face= 0; face< _faceNum; face++){
            float p = __mu[face][0];
            float r = __mu[face][1];
            float base = 1-p*p-r*r;
            float a_1 = (-(p-1)*(p-1)-r*r)/base;
            float a_2 = 2*r/base;
            float a_3 = a_2;   
            float a_4 = (-(1+p)*(1+p)-r*r)/base;
            // as they are in negative symmetry, the inverse is negative of the original  
            __muMatrix_U(face, face) = a_1;
            __muMatrix_U(face, _faceNum+face) = a_2;
            __muMatrix_U(_faceNum+face, face) = a_3;
            __muMatrix_U(_faceNum+face, _faceNum+face) = a_4;
            __muMatrix_V(face, face) = -a_1; 
            __muMatrix_V(face, _faceNum+face) = -a_2;
            __muMatrix_V(_faceNum+face, face) = -a_3;
            __muMatrix_V(_faceNum+face, _faceNum+face) = -a_4; 
        } 
        MatrixXf __fullGradientOnVertex_U = operators.fullGradientOnVertex;
        MatrixXf __fullGradientOnVertex_V = operators.fullGradientOnVertex;  
        // compose final least square matrix
        __lsqMatrix_U = operators.fullDivergenceOnFace* __muMatrix_U* __fullGradientOnVertex_U;
        __lsqMatrix_V = operators.fullDivergenceOnFace* __muMatrix_V* __fullGradientOnVertex_V; 
        // handle boundary conditions and organize least square vector
        VectorXf __lsqVector_U(__lsqMatrix_U.rows());
        VectorXf __lsqVector_V(__lsqMatrix_V.rows());
        __lsqVector_U.setZero();
        __lsqVector_V.setZero();  
        // horizontal boundary
        float b = __bounds[2];
        for (int &j : _vertexIndexAtRectBoundaryHorizontal_Top){  
            __lsqVector_V(j) = 0; 
            __lsqMatrix_V(j,j) = -1;
            __lsqVector_V += b* __lsqMatrix_V.col(j);
            for (int i=0; i<__lsqMatrix_V.rows(); i++){  
                    __lsqMatrix_V(i,j) = 0;   
            }for (int i=0; i<__lsqMatrix_V.cols(); i++){  
                    __lsqMatrix_V(j,i) = 0;    
            }__lsqMatrix_V(j,j) = 1; 
        }b = __bounds[3];
        for (int &j : _vertexIndexAtRectBoundaryHorizontal_Down){     
            __lsqVector_V(j) = 0; 
            __lsqMatrix_V(j,j) = -1;
            __lsqVector_V += b* __lsqMatrix_V.col(j);
            for (int i=0; i<__lsqMatrix_V.rows(); i++){  
                    __lsqMatrix_V(i,j) = 0;   
            }for (int i=0; i<__lsqMatrix_V.cols(); i++){  
                    __lsqMatrix_V(j,i) = 0;    
            }__lsqMatrix_V(j,j) = 1; 
        } 
        // vertical boundary
        b = __bounds[0];  
        for (int &j : _vertexIndexAtRectBoundaryVertical_Left){   
            __lsqVector_U(j) = 0;
            __lsqMatrix_U(j,j) = -1; 
            __lsqVector_U += b* __lsqMatrix_U.col(j); 
            for (int i=0; i<__lsqMatrix_U.rows(); i++){  
                    __lsqMatrix_U(i,j) = 0;   
            }for (int i=0; i<__lsqMatrix_U.cols(); i++){  
                    __lsqMatrix_U(j,i) = 0;    
            }__lsqMatrix_U(j,j) = 1; 
        }b = __bounds[1];
        for (int &j : _vertexIndexAtRectBoundaryVertical_Right){      
            __lsqVector_U(j) = 0;
            __lsqMatrix_U(j,j) = -1; 
            __lsqVector_U += b* __lsqMatrix_U.col(j);
            for (int i=0; i<__lsqMatrix_U.rows(); i++){  
                    __lsqMatrix_U(i,j) = 0;   
            }for (int i=0; i<__lsqMatrix_U.cols(); i++){  
                    __lsqMatrix_U(j,i) = 0;    
            }__lsqMatrix_U(j,j) = 1; 
        } 
        VectorXf sol_U = __lsqMatrix_U.jacobiSvd(ComputeThinU | ComputeThinV).solve(-__lsqVector_U);
        VectorXf sol_V = __lsqMatrix_V.jacobiSvd(ComputeThinU | ComputeThinV).solve(-__lsqVector_V); 
        // fill back solution
        for (int i = 0;i< _vertexNum;i++){
            __solution[i][0] = sol_U(i);
            __solution[i][1] = sol_V(i);
        }
    } 
} 
 
