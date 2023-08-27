///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  AttributedObject.cpp
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

// include the header file
#include "AttributedObject.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

// include the Cartesian 3- vector class
#include "Cartesian3.h"

#define MAXIMUM_LINE_LENGTH 1024
#define REMAP_TO_UNIT_INTERVAL(x) (0.5 + (0.5 * (x)))
#define REMAP_FROM_UNIT_INTERVAL(x) (-1.0 + (2.0 * (x)))

#define N_ITERATIONS 100000

// constructor will initialise to safe values
AttributedObject::AttributedObject()
    : centreOfGravity(0.0, 0.0, 0.0)
{ // AttributedObject()
    // force arrays to size 0
    vertices.resize(0);
    colours.resize(0);
    normals.resize(0);
    textureCoords.resize(0);
    firstDirectedEdge.resize(0);
    faceVertices.resize(0);
    otherHalf.resize(0);
    verts.resize(0);
    edges.resize(0);

} // AttributedObject()

// read routine returns true on success, failure otherwise
bool AttributedObject::ReadObjectStream(std::istream &geometryStream)
{ // ReadObjectStream()

    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];

    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
    { // not eof
        // character to read
        char firstChar = geometryStream.get();

        //         std::cout << "Read: " << firstChar << std::endl;

        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the character we read
        switch (firstChar)
        {         // switch on first character
        case '#': // comment line
            // read and discard the line
            geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            break;

        case 'v': // vertex data of some type
        {         // some sort of vertex data
            // retrieve another character
            char secondChar = geometryStream.get();

            // bail if we ran out of file
            if (geometryStream.eof())
                break;

            // now use the second character to choose branch
            switch (secondChar)
            {         // switch on second character
            case ' ': // space - indicates a vertex
            {         // vertex read
                Cartesian3 vertex;
                geometryStream >> vertex;
                vertices.push_back(vertex);
                //                         std::cout << "Vertex " << vertex << std::endl;
                break;
            }         // vertex read
            case 'c': // c indicates colour
            {         // normal read
                Cartesian3 colour;
                geometryStream >> colour;
                colours.push_back(colour);
                //                         std::cout << "Colour " << colour << std::endl;
                break;
            }         // normal read
            case 'n': // n indicates normal vector
            {         // normal read
                Cartesian3 normal;
                geometryStream >> normal;
                normals.push_back(normal);
                //                         std::cout << "Normal " << normal << std::endl;
                break;
            }         // normal read
            case 't': // t indicates texture coords
            {         // tex coord
                Cartesian3 texCoord;
                geometryStream >> texCoord;
                textureCoords.push_back(texCoord);
                //                         std::cout << "Tex Coords " << texCoord << std::endl;
                break;
            } // tex coord
            default:
                break;
            } // switch on second character
            break;
        } // some sort of vertex data

        case 'f': // face data
        {         // face
            // make a hard assumption that we have a single triangle per line
            unsigned int vertexID;

            // read in three vertices
            for (unsigned int vertex = 0; vertex < 3; vertex++)
            { // per vertex
                // read a vertex ID
                geometryStream >> vertexID;

                // subtract one and store them (OBJ uses 1-based numbering)
                faceVertices.push_back(vertexID - 1);
            } // per vertex
            break;
        } // face

            // default processing: do nothing
        default:
            break;

        } // switch on first character

    } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
    { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];

        // and K_divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();

            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;

        } // per vertex
    }     // non-empty vertex set

    // 	std::cout << "Centre of Gravity: " << centreOfGravity << std::endl;
    // 	std::cout << "Object Size:       " << objectSize << std::endl;

    // return a success code
    return true;
} // ReadObjectStream()

// write routine
void AttributedObject::WriteObjectStream(std::ostream &geometryStream)
{ // WriteObjectStream()
    geometryStream << "# " << (faceVertices.size() / 3) << " triangles" << std::endl;
    geometryStream << std::endl;

    // output the vertex coordinates
    geometryStream << "# " << vertices.size() << " vertices" << std::endl;
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "v  " << std::fixed << vertices[vertex] << std::endl;

    // output the vertex colours
    geometryStream << "# " << colours.size() << " vertex colours" << std::endl;
    for (unsigned int vertex = 0; vertex < colours.size(); vertex++)
        geometryStream << "vc " << std::fixed << colours[vertex] << std::endl;

    // output the vertex normals
    geometryStream << "# " << normals.size() << " vertex normals" << std::endl;
    for (unsigned int vertex = 0; vertex < normals.size(); vertex++)
        geometryStream << "vn " << std::fixed << normals[vertex] << std::endl;

    // output the vertex coords
    geometryStream << "# " << textureCoords.size() << " vertex tex coords" << std::endl;
    for (unsigned int vertex = 0; vertex < textureCoords.size(); vertex++)
        geometryStream << "vt " << std::fixed << textureCoords[vertex] << std::endl;

    // and the faces
    for (unsigned int face = 0; face < faceVertices.size(); face += 3)
    { // per face
        geometryStream << "f";

        // loop through # of vertices
        for (unsigned int vertex = 0; vertex < 3; vertex++)
        { // per vertex
            geometryStream << " ";
            geometryStream << faceVertices[face + vertex] + 1;
        } // per vertex
        // end the line
        geometryStream << std::endl;
    } // per face

} // WriteObjectStream()

// routine to render
void AttributedObject::Render(RenderParameters *renderParameters)
{ // Render()
    // make sure that textures are disabled
    glDisable(GL_TEXTURE_2D);

    float scale = renderParameters->zoomScale;
    
    scale /= objectSize;
    // Scale defaults to the zoom setting
    glTranslatef(-centreOfGravity.x * scale, -centreOfGravity.y * scale, -centreOfGravity.z * scale);

    if (renderParameters->useWireframe)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // start rendering
    glBegin(GL_TRIANGLES);

    // loop through the faces: note that they may not be triangles, which complicates life
    for (unsigned int face = 0; face < faceVertices.size(); face += 3)
    { // per face

        // now do a loop over three vertices
        for (unsigned int vertex = 0; vertex < 3; vertex++)
        { // per vertex
            // set colour using vertex ID
            // if (renderParameters->useTexCoords || renderParameters->useNormal)
            // {
            //     glColor3f(
            //                 textureCoords[faceVertices[face + vertex]].x * 255,
            //             textureCoords[faceVertices[face + vertex]].y * 255,
            //             textureCoords[faceVertices[face + vertex]].z * 255);
            // }
            // if (renderParameters->useNormal)
            // {
            //     glColor3f(
            //                 normals[faceVertices[face + vertex]].x * 255,
            //             normals[faceVertices[face + vertex]].y * 255,
            //             normals[faceVertices[face + vertex]].z * 255);
            // }
            // else
            // {
            // }

            // if (renderParameters->renderTexture || renderParameters->renderNormalMap)
            if (renderParameters->renderTexture && (!renderParameters->useTexCoords)) // Only Click the "Texture" CheckBox
            {
                glColor3f(
                        colours[faceVertices[face + vertex]].x,
                        colours[faceVertices[face + vertex]].y,
                        colours[faceVertices[face + vertex]].z);
                glVertex3f(
                        scale * textureCoords[faceVertices[face + vertex]].x,
                        scale * textureCoords[faceVertices[face + vertex]].y,
                        scale * textureCoords[faceVertices[face + vertex]].z);
            }
            else if (renderParameters->renderTexture && renderParameters->useTexCoords) // Click Both "Texture" CheckBox and "UVW -> RGB" CheckBox
            {
                glColor3f(
                        textureCoords[faceVertices[face + vertex]].x * 255,
                        textureCoords[faceVertices[face + vertex]].y * 255,
                        textureCoords[faceVertices[face + vertex]].z * 255);
                glVertex3f(
                        scale * textureCoords[faceVertices[face + vertex]].x,
                        scale * textureCoords[faceVertices[face + vertex]].y,
                        scale * textureCoords[faceVertices[face + vertex]].z);
            }
            else if (renderParameters->renderNormalMap && (!renderParameters->useNormal))
            {
                glColor3f(
                        normals[faceVertices[face + vertex]].x,
                        normals[faceVertices[face + vertex]].y,
                        normals[faceVertices[face + vertex]].z);
                glVertex3f(
                        scale * textureCoords[faceVertices[face + vertex]].x,
                        scale * textureCoords[faceVertices[face + vertex]].y,
                        scale * textureCoords[faceVertices[face + vertex]].z
                        );
            }
            else if (renderParameters->renderNormalMap && renderParameters->useNormal)
            {
                glColor3f(
                        normals[faceVertices[face + vertex]].x * 255,
                        normals[faceVertices[face + vertex]].y * 255,
                        normals[faceVertices[face + vertex]].z * 255);
                glVertex3f(
                        scale * textureCoords[faceVertices[face + vertex]].x,
                        scale * textureCoords[faceVertices[face + vertex]].y,
                        scale * textureCoords[faceVertices[face + vertex]].z
                        );
            }
            else // Without CheckBox Situation
            {
                glColor3f(
                        colours[faceVertices[face + vertex]].x,
                        colours[faceVertices[face + vertex]].y,
                        colours[faceVertices[face + vertex]].z);

                glVertex3f(
                        scale * vertices[faceVertices[face + vertex]].x,
                        scale * vertices[faceVertices[face + vertex]].y,
                        scale * vertices[faceVertices[face + vertex]].z);
            }
        } // per vertex
    }     // per face

    // close off the triangles
    glEnd();

    // revert render mode
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

// Enforce Half-Edge Data Structure
void AttributedObject::Pairing()
{
    // Init & Build Half Edges
    BuildHalfEdges();
    // Init & Build Opposite Edges
    BuildOtherHalfs();
    // Init & Build First Directed Edges
    BuildFirstDirectedEdges();
    // Init & Build Attributes of Vertices
    BuildAttributes();
}

// Total Process of Texture Parameterization Process, And it could be divided into 2 steps:
// 1. Fix Boundary as A Square (or a circle)
// 2. Using Iterative Method, Here using Averge Weight Mehthod to create a new uv image
void AttributedObject::TextureParameterization()
{
    FixBoundary();

    SetupInteralVertices();
}

// Find Neigbors of a Vertex
std::vector<unsigned int> AttributedObject::FindNeigbors(unsigned int ID)
{
    // get all vector of adjacent vertex
    std::vector<unsigned int> neigbors;

    // get all vector of adjacent vertex
    const int &start_edge = verts[ID].fd_edge;
    int curr_edge = start_edge;
    if (start_edge == -1)
    {
        return neigbors;
    }
    do
    {
        neigbors.emplace_back(edges[PREV(curr_edge)].vertID);
        curr_edge = edges[curr_edge].oppo;
        if (edges[PREV(curr_edge)].vertID == ID)
        {
            curr_edge = edges[curr_edge].next;
        }
        else
        {
            break;
        }
    } while (curr_edge != start_edge);

    return neigbors;
}

void AttributedObject::FixBoundary()
{
    // Find Vertices of a Hole at Speicfic Sequence
    std::vector<unsigned int> boundaryVerts = FindHoleVerts();

    // From the Boundaries Fix them into a A Square
    FixToASquare(boundaryVerts);
}

// Find Vertices of a Hole at Speicfic Sequence
std::vector<unsigned int> AttributedObject::FindHoleVerts()
{

    // Convert to conected boundary vertices
    std::vector<unsigned int> unique_boundaryVerts;

    const auto &start_vert = edges[boundaryEdges[0]].vertID; // start vert of the search loop, if combines as a loop circle, immediately jump out.

    unique_boundaryVerts.emplace_back(start_vert);

    bool stop = false;                         // flag of jump out of the loop
    unsigned int count = 1;                    // count of searched boundary edges
    unsigned int curr_edge = boundaryEdges[0]; // edge in searching, initalize it with index 0

    while (count < boundaryEdges.size() && stop == false)
    {
        for (unsigned int j = 0; j < boundaryEdges.size() && stop == false; j++)
        {
            const auto &search_vert = edges[PREV(curr_edge)].vertID; // Using aliasing to make lifies easier
            const auto &search_edge = edges[boundaryEdges[j]];

            // I Expect the searched edge's start ID = current edge Point ID, So they could Link Togther.
            if (search_vert == search_edge.vertID && search_vert == start_vert)
            {
                verts[start_vert].boundary = true;
                verts[search_vert].boundary = true;

                stop = true;
            }
            else if (search_vert == search_edge.vertID && search_vert != start_vert)
            {
                // Acculuate count number as the flag of jumping out the loop.
                count++;
                curr_edge = boundaryEdges[j];

                unique_boundaryVerts.emplace_back(edges[boundaryEdges[j]].vertID);
                verts[search_edge.vertID].boundary = true;

                break;
            }
        }
    }

    return unique_boundaryVerts;
}

// From the Boundaries Fix them into a A Square
void AttributedObject::FixToASquare(std::vector<unsigned int> boundaries)
{
    float K_sum;
    for (unsigned int v = 0; v < boundaries.size() - 1; ++v)
        K_sum += (verts[boundaries[v + 1]].position - verts[boundaries[v]].position).length();

    for (unsigned int u = 0; u < boundaries.size(); u++)
    {
        float K_sub = 0;

        if (u >= 1)
            for (unsigned int w = 0; w < u - 1; ++w)
                K_sub += (verts[boundaries[w + 1]].position - verts[boundaries[w]].position).length();

        float K_div = K_sub / K_sum;

        if (K_div >= 0 && K_div < 0.25)
            textureCoords[boundaries[u]] = vec3(0.5 - 4 * K_div, 1, 0);
        else if (K_div >= 0.25 && K_div < 0.5)
            textureCoords[boundaries[u]] = vec3(-0.5, 2 - 4 * K_div, 0);
        else if (K_div >= 0.5 && K_div < 0.75)
            textureCoords[boundaries[u]] = vec3(4 * K_div - 2.5, 0, 0);
        else if (K_div >= 0.75 && K_div < 1)
            textureCoords[boundaries[u]] = vec3(0.5, 4 * K_div - 3, 0);
    }
}

void AttributedObject::SetupInteralVertices()
{
    // Init Interior Vertices to the Center (from Slides)
    for (unsigned int i = 0; i < verts.size(); i++)
        if (verts[i].boundary == false)
            textureCoords[i] = vec3(0.0f, 0.5f, 0.0f);

    // Iterative to Get Results
    for (unsigned int depth = 0; depth <= 150; ++depth)
        // Averge Weight Mehthod to create a new uv image
        AvgWeightMethod();
}

// Averge Weight Mehthod to create a new uv image
void AttributedObject::AvgWeightMethod()
{
    for (unsigned int ID = 0; ID < verts.size(); ID++)
    {// Iterative Verts list and Compute the New position depends on neigbor weights.  
        if (verts[ID].boundary == false)
        {
            std::vector<unsigned int> neigbors = FindNeigbors(ID);
            if (neigbors.size() == 0)
            {
                textureCoords[ID] = vec3(0.0f, 0.5f, 0.0f);
            }
            vec3 nc_sum;
            for (const auto &n : neigbors)
            {
                nc_sum += textureCoords[n];
            }
            textureCoords[ID] =  nc_sum / neigbors.size();
        }
    }
}

// Update Smoothed Normals and Assign Colors from UV
void AttributedObject::UpdateNormalsFromUV()
{
    // Clear Normals from exceptions
    this->normals.clear();
    for (unsigned int v = 0; v < verts.size(); ++v)
    {
        // GetNeigbors list
        std::vector<unsigned int> neigbors = FindNeigbors(v);

        // Smooth Weights Parameters
        float sum_area = 0.0f;
        vec3 sum_normal = vec3();
        
        // Exclude with Neighbors size < 2 ---> That means Boundaries. 
        if (neigbors.size() > 2)
        {
            for (unsigned int n = 0; n < neigbors.size() - 1; ++n)
            {
                // other two vert IDs
                const auto& start = neigbors[n];
                const auto& point = neigbors[n + 1];

                // Three Positions
                const auto& v0 = verts[v].position;
                const auto& v1 = verts[start].position;
                const auto& v2 = verts[point].position;

                // edges of vertices
                const auto e0 = v0 - v1;
                const auto e1 = v1 - v2;

                vec3 face_normal = (e0).cross(e1).unit();
                sum_normal += face_normal;
                float area = (e0).cross(e1).dot(face_normal);
                sum_area += area;
            }
            // auto normal = sum_normal / sum_area;
            this->normals.emplace_back(sum_normal / sum_area);
        }
        else
            this->normals.emplace_back(vec3(0, 0, 1)); // just return no color
    }
}

// Init & Build Half Edges
void AttributedObject::BuildHalfEdges()
{
    for (auto f = 0; f < this->faceVertices.size() / 3; ++f)
    {
        //                     vertID      next       prev       face
        HEEdge e0(faceVertices[f * 3    ], f * 3 + 1, f * 3 + 2, f);
        HEEdge e1(faceVertices[f * 3 + 1], f * 3 + 2, f * 3,     f);
        HEEdge e2(faceVertices[f * 3 + 2], f * 3    , f * 3 + 1, f);

        edges.emplace_back(e0);
        edges.emplace_back(e1);
        edges.emplace_back(e2);
    }
}

// Init & Build Opposite Edges
void AttributedObject::BuildOtherHalfs()
{
    for (unsigned int index = 0; index < edges.size(); ++index)
    {
        // v_start
        const auto &v_start = edges[index].vertID;

        // v_point
        const auto &v_point = edges[PREV(index)].vertID;

        bool isBorder = true;
        for (unsigned int inner = 0; inner < edges.size(); ++inner)
        {
            const unsigned int &i_start = edges[inner].vertID;
            const unsigned int &i_point = edges[PREV(inner)].vertID;
            if (i_start == v_point && i_point == v_start)
            {
                edges[index].oppo = inner;
                isBorder = false;
                break;
            }
        }

        if (isBorder)
        {
            boundaryEdges.emplace_back(index);
            edges[index].edge_boundary = true;
        }
    }
}

// Init & Build First Directed Edges
void AttributedObject::BuildFirstDirectedEdges()
{
    // Update First Directed Edges
    for (auto vertID = 0; vertID < verts.size(); ++vertID)
    {
        bool stop = false;
        for (auto edgeID = 0; edgeID < edges.size(); ++edgeID)
        {
            auto &start_vert = edges[edgeID].vertID;
            if (start_vert == vertID && stop == false)
            {
                verts[vertID].fd_edge = edgeID;
                stop = true;
            }
        }
    }
}

// Init & Build Attributes of Vertices
void AttributedObject::BuildAttributes()
{
    for (unsigned int v = 0; v < vertices.size(); v++)
    {
        verts.emplace_back(vertices[v]);
        textureCoords.emplace_back(vec3(0.0f, 0.5f, 0.0f));
    }

    // Compute Face Normals Here, ignore here
    // for (unsigned int f = 0; f < faceVertices.size() / 3; ++f)
    // {
    //     Cartesian3 e1 = vertices[3 * f + 2] - vertices[3 * f];
    //     Cartesian3 e2 = vertices[3 * f + 1] - vertices[3 * f];
    //     Cartesian3 normal = e1.cross(e2).unit();
    //     normals.emplace_back(normal);
    // }
}
