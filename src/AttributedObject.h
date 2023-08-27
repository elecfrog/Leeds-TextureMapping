///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  AttributedObject.h
//  ------------------------
//
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and
//  rendering an object file
//
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or
//  normal direction, &c.  And if it doesn't work on
//  all object files, that's fine.
//
//	Variant on TexturedObject that stores explicit RGB
//	values for each vertex
//
///////////////////////////////////////////////////

// include guard for AttributedObject
#ifndef _ATTRIBUTED_OBJECT_H
#define _ATTRIBUTED_OBJECT_H

// include the C++ standard libraries we need for the header
#include <vector>
#include <iostream>
#include <array>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

// include the unit with Cartesian 3-vectors
#include "Cartesian3.h"
// the render parameters
#include "RenderParameters.h"

// define a macro for "not used" flag
#define NO_SUCH_ELEMENT -1

// use macros for the "previous" and "next" IDs
#define NEXT(x) ((x) % 3) ? ((x) - 1) : ((x) + 2)
#define PREV(x) (((x) % 3) == 2) ? ((x) - 2) : ((x) + 1)

// Vertex Definition
struct HEVertex {

    vec3 position;          // Vertex position
    int fd_edge;            // First directed edge
    bool boundary;          // The flag 

    HEVertex() = default; // Set as default constructor
    
    HEVertex(vec3 _position) // Parameter constructor
        :position(_position), fd_edge(-1), boundary(false) { }
};

// Half Edge Definition
struct HEEdge {

    int vertID;         // start from vert ID
    int next;           // next edge IDs
    int prev;           // prev edge IDs
    int oppo;           // oppo edge IDs
    int face;           // face IDs
    bool edge_boundary; // Flag of it is a boundary

    HEEdge() = default; // Set as default constructor
    
    HEEdge(int _vertID, int _next, int _prev, int _face) // Parameter constructor
        :vertID(_vertID), next(_next), prev(_prev), face(_face), oppo(-1), edge_boundary(false) { }
};

class AttributedObject
    { // class AttributedObject
    public:
    // vector of vertices
    std::vector<Cartesian3> vertices;

    // vector of colours stored as cartesian triples in float
    std::vector<Cartesian3> colours;

    // vector of normals
    std::vector<Cartesian3> normals;

    // vector of normals of faces
    std::vector<Cartesian3> face_normals;

    // vector of texture coordinates (stored as triple to simplify code)
    std::vector<Cartesian3> textureCoords;

    // vector for the "first" directed edge for each vertex
    std::vector<long> firstDirectedEdge;

    // vector of faces - doubles as the "to" array for edges
    std::vector<unsigned int> faceVertices;

    // vector of the "other half" of the directed edges
    std::vector<long> otherHalf;

    // centre of gravity - computed after reading
    Cartesian3 centreOfGravity;

    // size of object - i.e. radius of circumscribing sphere centred at centre of gravity
    float objectSize;

    // Vertices (with Attributs) Containter
    std::vector<HEVertex> verts;

    // Half Edges Containter
    std::vector<HEEdge> edges;

    // constructor will initialise to safe values
    AttributedObject();

    // read routine returns true on success, failure otherwise
    bool ReadObjectStream(std::istream &geometryStream);

    // write routine
    void WriteObjectStream(std::ostream &geometryStream);

    // routine to render
    void Render(RenderParameters *renderParameters);

/***************************************
* TASK I: Build HalfEdge Data Structure
****************************************/
private:
    // Init & Build Half Edges
    void BuildHalfEdges();

    // Indices of Edge Boundaries
    std::vector<unsigned int> boundaryEdges;
    
    // Init & Build Opposite Edges
    void BuildOtherHalfs();

    // Init & Build First Directed Edges
    void BuildFirstDirectedEdges();
    
    // Init & Build Attributes of Vertices
    void BuildAttributes();
public:
    // Enforce Half-Edge Data Structure
    void Pairing();

/***************************************
* TASK II: Identify Boundaries for Floater's Algorithm
****************************************/
private:
    // Find Vertices of a Hole at Speicfic Sequence
    std::vector<unsigned int> FindHoleVerts();
    // From the Boundaries Fix them into a A Square
    void FixToASquare(std::vector<unsigned int> _boundaries);
    // Entry of Boundary Operations
    void FixBoundary();


/*****************************************************
* TASK III: TextureParameterization using Floater's Algorithm
******************************************************/
private:
    void SetupInteralVertices();
    // Averge Weight Mehthod to create a new uv image
    void AvgWeightMethod();
    // Find Neigbors of a Vertex
    std::vector<unsigned int> FindNeigbors(unsigned int ID);
public:
    // Entry of Texture Parameterization Process, And it could be divided into 2 steps:
    // 1. Fix Boundary as A Square (or a circle)
    // 2. Using Iterative Method, Here using Averge Weight Mehthod to create a new uv image
    void TextureParameterization();

/*****************************************************
* TASK VII: Normal Map
******************************************************/
    // Update Smoothed Normals and Assign Colors from UV
    void UpdateNormalsFromUV();


}; // class AttributedObject

// end of include guard for AttributedObject
#endif
