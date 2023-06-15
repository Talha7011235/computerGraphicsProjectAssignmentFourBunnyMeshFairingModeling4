// Talha Hussain
// Computer Graphics Project Assignment 4 Mesh Fairing Modeling
// main4.cpp
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "io.h"
#include "LinearSolver/cghs.h"

typedef std::vector<Eigen::Vector3d> vertex_list;
typedef std::vector<Tri> triangle_list;

// The Structure Data Type "struct Operator" is used to perform linear solving using the cghs.h Header File from the LinearSolver package directory folder.
struct Operator
{
    int n;
    double lambda, dt, eps, strength;
    bool biharmonic;
    triangle_list triangles;
};

// The Strucutre Data Type "struct LoopSubdivision" keeps track of the newly subdivided model.
struct LoopSubdivision {
    vertex_list vertices;
    std::vector<Tri> triangles;
};

// ****************************************************************************************************************************************************************
// The Definiton of Function "bool is_digits(const std::string &str)" is that the Function "bool is_digits(const std::string &str)" takes in a string and returns *
// whether it consists of only digits.                                                                                                                            *
// ****************************************************************************************************************************************************************
bool is_digits(const std::string &str)
{
    return str.find_first_not_of("0123456789.") == std::string::npos;
}

// The Function "vertex_list basicLaplacianSmoothing (const vertex_lists &points, const triangle_list &triangles, const double &strength, const bool &biharmonic)"
// performs basic Laplacian Smoothing as a result of using the umbrella or biharmonic operator.
vertex_list basicLaplacianSmoothing(
    const vertex_list &points,
    const triangle_list &triangles,
    const double &strength,
    const bool &biharmonic)
{
    size_t numVertices = points.size();
    vertex_list verticies;
    std::vector<size_t> edgeSet;
    verticies.resize(numVertices, Eigen::Vector3d::Zero());
    edgeSet.resize(numVertices, 0);
    
    for (size_t e = 0; e < triangles.size(); ++e)
    {
        for(size_t j = 0; j < 3; ++j) {
            // Find the point.
            int vertex = triangles[e].indices[j];
            verticies[vertex] += (points[triangles[e].indices[(j + 1) % 3]] - points[vertex]);
            verticies[vertex] += (points[triangles[e].indices[(j + 2) % 3]] - points[vertex]);
            edgeSet[vertex] += 2;
        }
    }

    for (size_t i = 0; i < numVertices; ++i)
    {
        Eigen::Vector3d pos = verticies[i];
        size_t numEdges = edgeSet[i];

        const double weight = pow(1.0 / numEdges, 2);
        if (biharmonic)
        {
            verticies[i] = points[i] + pow(strength * weight, 2.0) * pos;
        }
        else
        {
            verticies[i] = points[i] + (strength * weight * pos);
        }
    }

    return verticies;
}

// The Function std::vector<std:map<size_t, double>> getCotangentWeigths(std::vector<std::map<size_t, double>> getCotangentWeights(const vertex_list &verticies,
// const triangle_list &triangles) calculates and returns a Vector of Cotangent Weights for use in Laplacian Smoothing.
std::vector<std::map<size_t, double>> getCotanWeights(
    const vertex_list &vertices,
    const triangle_list &triangles)
{
    size_t numVertices = vertices.size();
    std::vector<std::map<size_t, double>> weights(numVertices);
    weights.resize(numVertices);

    size_t triCount = triangles.size();

    for (size_t e = 0; e < triCount; ++e)
    {
        for(int j = 0; j < 3; ++j) {
            // Find the point.
            int index = triangles[e].indices[j];
            int next = triangles[e].indices[(j + 1) % 3];
            int prev = triangles[e].indices[(j + 2) % 3];            

            Eigen::Vector3d A = vertices[index] - vertices[prev];
            Eigen::Vector3d B = vertices[next] - vertices[prev];
            Eigen::Vector3d C = vertices[index] - vertices[next];
            Eigen::Vector3d D = vertices[prev] - vertices[next];

            double cotangentAB = A.dot(B) / (std::numeric_limits<double>::epsilon() + A.cross(B).norm());
            double cotangentCD = C.dot(D) / (std::numeric_limits<double>::epsilon() + C.cross(D).norm());

            if(weights[index].find(next) == weights[index].end()) {
                weights[index][next] = 0.0;
            }
            
            if(weights[index].find(prev) == weights[index].end()) {
                weights[index][prev] = 0.0;
            }

            weights[index][next] += cotangentAB * 0.5;
            weights[index][prev] += cotangentCD * 0.5;
        }
    }

    return weights;
}

// The Function "vertex_list laplacianSmoothingWithCotangents(const vertex_list &points, const triangle_list &triangles, const size_t &strength, const
// bool &biharmonic)" performs Laplacian Smoothing using the Cotangent Weights instead of using the simple accumulation.
vertex_list laplacianSmoothingWithCotangents(
    const vertex_list &points,
    const triangle_list &triangles,
    const size_t &strength,
    const bool &biharmonic)
{
    size_t numVertices = points.size();
    size_t numTriangles = triangles.size();
    vertex_list verticies;
    std::vector<double> sum;
    verticies.resize(numVertices, Eigen::Vector3d::Zero());
    sum.resize(numVertices, 0.0);
    
    std::vector<std::map<size_t, double>> weights = getCotanWeights(points, triangles);
    
    for (size_t e = 0; e < numTriangles; ++e)
    {
        for(int j = 0; j < 3; ++j) {
            // Find the point.
            int index = triangles[e].indices[j];
            int next = triangles[e].indices[(j + 1) % 3];
            int prev = triangles[e].indices[(j + 2) % 3];

            double nextWeight = weights[index][next];
            double prevWeight = weights[index][prev];
            if (biharmonic)
            {
                nextWeight = pow(nextWeight, 2.0);
                prevWeight = pow(prevWeight, 2.0);
            }

            verticies[index] += points[next] * (nextWeight * 0.5);
            verticies[index] += points[prev] * (prevWeight * 0.5);
            sum[index] += (nextWeight + prevWeight) * 0.5;
        }
    }

    for (size_t i = 0; i < numVertices; ++i)
    {
        verticies[i] = (verticies[i] / sum[i]) * strength + points[i] * (1.0-strength);
    }

    return verticies;
}

// ************************************************************************************************************************************
// The Definition of Function "double getInternalVolume(const vertex_list &verticies, const std::vector<Tri> &triangles)" is that the *
// Function "double getInternalVolume(const vertex_list &verticies, const std::vector<Tri> &triangles" calculates the Internal Volume *
// of the object used for Volume Preservation.                                                                                        *
// ************************************************************************************************************************************
double getInternalVolume(const vertex_list &vertices, const std::vector<Tri> &triangles)
{
    size_t numTriangles = triangles.size();
    double volume = 0;
    for (size_t i = 0; i < numTriangles; ++i)
    {
        Tri triangle = triangles[i];
        Eigen::Vector3d center = (vertices[triangle.indices[0]] + vertices[triangle.indices[1]] + vertices[triangle.indices[2]]) / 3;
        Eigen::Vector3d normal = (vertices[triangle.indices[0]] - vertices[triangle.indices[1]]).cross(vertices[triangle.indices[0]] - vertices[triangle.indices[2]]);

        volume += center.dot(normal);
    }

    return volume / 6.0;
}

// ****************************************************************************************************
// The Definition of Function "void scaleVerticies(vertex_list &verticies, double scale)" is that the *
// Function "void scaleVerticies(vertex_list &verticies, double scale)" uniformly scales the object,  *
// used to preserve the Internal Volume of the object.                                                *
// ****************************************************************************************************
void scaleVertices(vertex_list &vertices, double scale)
{
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        vertices[i] *= scale;
    }
}

// ***********************************************************************************************************************************************************
// The Definiton of Function "void mult(const Operator &op, double *v, double *w)" is that the Function void mult(const Operator &op, double *v, double *w)" *
// performs the Linear Solving Multiplication Operation using the Operation object and list of input and output.                                             *
// ***********************************************************************************************************************************************************
void mult(const Operator &op, double *v, double *w)
{
    vertex_list pts;
    pts.resize(op.n);
    for (unsigned int i = 0; i < pts.size(); i++)
    {
        pts[i][0] = v[3 * i + 0];
        pts[i][1] = v[3 * i + 1];
        pts[i][2] = v[3 * i + 2];
    }
    vertex_list l = basicLaplacianSmoothing(pts, op.triangles, op.strength, op.biharmonic);
    for (unsigned int i = 0; i < pts.size(); i++)
    {
        l[i] *= op.lambda * op.dt;
        w[3 * i + 0] = v[3 * i + 0] - l[i][0];
        w[3 * i + 1] = v[3 * i + 1] - l[i][1];
        w[3 * i + 2] = v[3 * i + 2] - l[i][2];
    }
}

// The Function "vertex_list performImplicit(const vertex_list &points, const Operator &op)" performs the Implicit Operation on the points using the
// Linear Solver Operation.
vertex_list performImplicit(const vertex_list &points, const Operator &op)
{
    size_t pointCount = points.size();
    double b[3 * pointCount];
    double x[3 * pointCount];

    for (size_t i = 0; i < pointCount; ++i)
    {
        b[3 * i + 0] = points[i][0];
        b[3 * i + 1] = points[i][1];
        b[3 * i + 2] = points[i][2];
        x[3 * i + 0] = points[i][0];
        x[3 * i + 1] = points[i][1];
        x[3 * i + 2] = points[i][2];
    }

    cghs<Operator>(3 * pointCount, op, b, x, op.eps, false);

    vertex_list output;
    output.resize(pointCount, Eigen::Vector3d::Zero());

    for (unsigned int i = 0; i < pointCount; i++)
    {
        output[i][0] = x[3 * i + 0];
        output[i][1] = x[3 * i + 1];
        output[i][2] = x[3 * i + 2];
    }
    return output;
}

// Implement BlasLib functions to allow for proper compiliation.

// **************************************************************************************************************************
// The Definition of Function "void dscal_(const int *n, const double *alpha, double *x, const int *incx" is that the       *
// Function "void dscal_(const int *n, const double *alpha, double *x, const int *incx" scales the n doubles in the x array *
// by the alpha value.                                                                                                      *
// **************************************************************************************************************************
void dscal_(const int *n, const double *alpha, double *x, const int *incx)
{
    // Linear scaling of a vector
    for (int i = 0; i < *n; ++i)
    {
        x[i * (*incx)] *= *alpha;
    }
}

// ****************************************************************************************************************************************
// The Definition of Function "void dcopy_(const int *n, const double *x, const int *incx, double *y, const int *incy)" is that the       *
// Function "void dcopy_(const int *n, const double *x< const int *incx, double *y, const int *incy)" copies the n doubles in the x array *
// into the y array along incremental values.                                                                                             *
// ****************************************************************************************************************************************
void dcopy_(const int *n, const double *x, const int *incx, double *y, const int *incy)
{
    // Copy Vector X onto vector Y
    for (int i = 0; i < *n; ++i)
    {
        y[i * (*incy)] = x[i * (*incx)];
    }
}

// *********************************************************************************************************************************************************
// The Definition of Function "void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *Y, const int *incy)" is that the   *
// Function "void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy)" performs a scale and translates *
// on the x array and then stores in the y array.                                                                                                          *
// *********************************************************************************************************************************************************
void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy)
{
    // y = a*x + y
    for (int i = 0; i < *n; ++i)
    {
        y[i * (*incy)] = *alpha * x[i * (*incx)] + y[i * (*incy)];
    }
}

// *****************************************************************************************************************************************
// The Definition of Function "double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy)" is that the *
// Function "double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy) performs a Dot Product on all  *
// of the values in the x and y arrays.                                                                                                    *
// *****************************************************************************************************************************************
double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy)
{
    double sum = 0;
    for (int i = 0; i < *n; ++i)
    {
        sum += x[i * (*incx)] * y[i * (*incy)];
    }
    return sum;
}
/** Finished with implementing the BlasLib functions to allow for proper complilation.*/

// The Function "LoopSubdivision Loopdivide(const vertex_list& vertices, const std::vector<Tri>& triangles) performs the Loop Subdivision operation and
// returns a new mesh object that has been subdivided and ready for smoothing.
LoopSubdivision Loopdivide(const vertex_list& vertices, const std::vector<Tri>& triangles) {
    LoopSubdivision output;
    std::vector<std::map<int, size_t>> edgeVerts;

    output.vertices.reserve(vertices.size() * 2);
    edgeVerts.resize(vertices.size());
    std::copy(vertices.begin(), vertices.end(), std::back_inserter(output.vertices));

    output.triangles.reserve(triangles.size() * 4);

    size_t triangleCount = triangles.size();
    for(size_t i = 0; i < triangleCount; ++i) {
        std::cout << "Generating Divisions: " << (i + 1) << "/" << triangleCount << std::endl
                  << "\x1b[A";
        Tri triangle = triangles[i];
        // Generate a center point for each edge of the Triangle.
        size_t center[] = { 0, 0, 0 };
        
        for(size_t j = 0; j < 3; ++j) {
            int min = std::min(triangle.indices[j], triangle.indices[(j+1) % 3]);
            int max = std::max(triangle.indices[j], triangle.indices[(j+1) % 3]);
            
            std::map<int, size_t> *edgeVert = &edgeVerts[min];
            if(edgeVert->find(max) != edgeVert->end()) {
                // If it does not exist, then add it to the set. 
                center[j] = (*edgeVert)[max];
                continue;
            }

            // It does not exist. Therefore, add it to the system.
            Eigen::Vector3d point = vertices[min] * 0.5 + vertices[max] * 0.5;
            center[j] = output.vertices.size();
            output.vertices.push_back(point);
            (*edgeVert)[max] = center[j];
        }

        /*    0
             / \
            c2-c0
           / \ / \
          2-- c1--1 */

        output.triangles.push_back(Tri(center[0], center[1], center[2]));
        output.triangles.push_back(Tri(triangle[0], center[0], center[2]));
        output.triangles.push_back(Tri(center[0], triangle[1], center[1]));
        output.triangles.push_back(Tri(center[2], center[1], triangle[2]));

        // For the following block of C++ code shown above, connect center points into new triangle and connect center points with original points to make
        // three triangles.
    }

    return output;
}

int main(int argc, char **argv)
{
    if (argc < 5)
    {
        std::cout << "USAGE:" << std::endl;
        std::cout << argv[0] << " [-c] input output stepsize number_of_iterations" << std::endl;
    }

    // Do Argument Parsing.
    bool useCotangent = false;
    bool preserveVolume = false;
    bool useBiHarmonic = false;
    bool useImplicit = false;
    bool useSubdivision = false;
    double eps = std::numeric_limits<double>::epsilon();
    int subdivisions = 0;
    int stage = 0;

    // Default Input Values.
    std::string strings[] = {
        "bunny.obj",
        "output.obj",
        "1",
        "1"};

    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "-b") == 0)
        {
            useBiHarmonic = true;
            continue;
        }
        if (strcmp(argv[i], "-c") == 0)
        {
            useCotangent = true;
            continue;
        }
        if (strcmp(argv[i], "-i") == 0)
        {
            useImplicit = true;
            i++;
            eps = atof(argv[i]);
            continue;
        }
        if (strcmp(argv[i], "-s") == 0)
        {
            useSubdivision = true;
            i++;
            subdivisions = atoi(argv[i]);
            continue;
        }
        if (strcmp(argv[i], "-v") == 0)
        {
            preserveVolume = true;
            continue;
        }
        if (stage < 4)
        {
            strings[stage++] = argv[i];
        }
    }

    // Ensure that the strength and number of iterations is correct.
    if (!is_digits(strings[2]))
    {
        std::cerr << "Invalid strength value" << std::endl;
        return 1;
    }

    if (!is_digits(strings[3]))
    {
        std::cerr << "Invalid number of iterations" << std::endl;
        return 1;
    }

    double strength = std::stof(strings[2]);
    size_t numIterations = std::stoul(strings[3]);

    vertex_list vertices;
    triangle_list triangles;

    // Read the input file into the Verticies List that is the "vertex_list" and Triangles List that is the "triangle_list".
    if (!readObjFile(const_cast<char *>(strings[0].c_str()), vertices, triangles))
    {
        std::cerr << "Unable to read object file" << std::endl;
        return -1;
    }

    // Perform any requested Subdivisions before processing list.
    if(useSubdivision) {
        for(int i = 0; i < subdivisions; ++i) {
            LoopSubdivision division = Loopdivide(vertices, triangles);
            vertices.clear();
            triangles.clear();
            vertices.reserve(division.vertices.size());
            triangles.reserve(division.triangles.size());
            std::copy(division.vertices.begin(), division.vertices.end(), std::back_inserter(vertices));
            std::copy(division.triangles.begin(), division.triangles.end(), std::back_inserter(triangles));
        }
        std::cout << std::endl;
    }

    // Retrieve the Internal Volume in order to use it for Volume Preservation.
    double volumeStart = getInternalVolume(vertices, triangles);
    double volumeEnd;
    Operator op;
    
    // If an implicit solver is desired, set up the operator.
    if(useImplicit) {
        op.n = vertices.size();
        op.strength = strength;
        op.lambda = 0.9;
        op.dt = 0.9;
        op.eps = eps;
        op.biharmonic = useBiHarmonic;
        op.triangles = triangles;
    }

    // Perform the smoothing for the requested number of iterations.
    for (size_t i = 0; i < numIterations; ++i)
    {
        std::cout << "Running Iterations: " << (i + 1) << "/" << numIterations << std::endl
                  << "\x1b[A";

        // Perform the desired smoothing operation.
        if (useImplicit)
        {
            vertices = performImplicit(vertices, op);
        }
        else if (useCotangent)
        {
            vertices = laplacianSmoothingWithCotangents(vertices, triangles, strength, useBiHarmonic);
        }
        else
        {
            vertices = basicLaplacianSmoothing(vertices, triangles, strength, useBiHarmonic);
        }

        // If Volume Preservation is desired, calculate and scale the object properly.
        if (preserveVolume)
        {
            volumeEnd = getInternalVolume(vertices, triangles);
            scaleVertices(vertices, pow(volumeStart / volumeEnd, 1 / 3.0));
        }
    }
    std::cout << std::endl;

    // Write the final object to the appropriate output file.
    std::cout << "Saving Output" << std::endl;
    writeObjFile(const_cast<char *>(strings[1].c_str()), vertices, triangles);

    return 0;
}
