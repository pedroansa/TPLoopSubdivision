#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <vector>
#include <memory>
#include <algorithm>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include <map>
#include <set>
#include <iostream>

class Mesh {
public:
    virtual ~Mesh();

    const std::vector<glm::vec3>& vertexPositions() const { return _vertexPositions; }
    std::vector<glm::vec3>& vertexPositions() { return _vertexPositions; }

    const std::vector<glm::vec3>& vertexNormals() const { return _vertexNormals; }
    std::vector<glm::vec3>& vertexNormals() { return _vertexNormals; }

    const std::vector<glm::vec2>& vertexTexCoords() const { return _vertexTexCoords; }
    std::vector<glm::vec2>& vertexTexCoords() { return _vertexTexCoords; }

    const std::vector<glm::uvec3>& triangleIndices() const { return _triangleIndices; }
    std::vector<glm::uvec3>& triangleIndices() { return _triangleIndices; }

    /// Compute the parameters of a sphere which bounds the mesh
    void computeBoundingSphere(glm::vec3& center, float& radius) const;

    void recomputePerVertexNormals(bool angleBased = false);
    void recomputePerVertexTextureCoordinates();

    void init();
    void initOldGL();
    void render();
    void clear();

    void addPlan(float square_half_side = 1.0f);

    void subdivideLinear() {
        std::vector<glm::vec3> newVertices = _vertexPositions;
        std::vector<glm::uvec3> newTriangles = _triangleIndices;

        struct Edge {
            unsigned int a, b;
            Edge(unsigned int c, unsigned int d) : a(std::min<unsigned int>(c, d)), b(std::max<unsigned int>(c, d)) {}
            bool operator < (Edge const& o) const { return a < o.a || (a == o.a && b < o.b); }
            bool operator == (Edge const& o) const { return a == o.a && b == o.b; }
        };
        std::map< Edge, unsigned int > newVertexOnEdge;
        for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
            unsigned int a = _triangleIndices[tIt][0];
            unsigned int b = _triangleIndices[tIt][1];
            unsigned int c = _triangleIndices[tIt][2];


            Edge Eab(a, b);
            unsigned int oddVertexOnEdgeEab = 0;
            if (newVertexOnEdge.find(Eab) == newVertexOnEdge.end()) {
                newVertices.push_back((_vertexPositions[a] + _vertexPositions[b]) / 2.f);
                oddVertexOnEdgeEab = newVertices.size() - 1;
                newVertexOnEdge[Eab] = oddVertexOnEdgeEab;
            }
            else { oddVertexOnEdgeEab = newVertexOnEdge[Eab]; }


            Edge Ebc(b, c);
            unsigned int oddVertexOnEdgeEbc = 0;
            if (newVertexOnEdge.find(Ebc) == newVertexOnEdge.end()) {
                newVertices.push_back((_vertexPositions[b] + _vertexPositions[c]) / 2.f);
                oddVertexOnEdgeEbc = newVertices.size() - 1;
                newVertexOnEdge[Ebc] = oddVertexOnEdgeEbc;
            }
            else { oddVertexOnEdgeEbc = newVertexOnEdge[Ebc]; }

            Edge Eca(c, a);
            unsigned int oddVertexOnEdgeEca = 0;
            if (newVertexOnEdge.find(Eca) == newVertexOnEdge.end()) {
                newVertices.push_back((_vertexPositions[c] + _vertexPositions[a]) / 2.f);
                oddVertexOnEdgeEca = newVertices.size() - 1;
                newVertexOnEdge[Eca] = oddVertexOnEdgeEca;
            }
            else { oddVertexOnEdgeEca = newVertexOnEdge[Eca]; }

            // set new triangles :
            newTriangles.push_back(glm::uvec3(a, oddVertexOnEdgeEab, oddVertexOnEdgeEca));
            newTriangles.push_back(glm::uvec3(oddVertexOnEdgeEab, b, oddVertexOnEdgeEbc));
            newTriangles.push_back(glm::uvec3(oddVertexOnEdgeEca, oddVertexOnEdgeEbc, c));
            newTriangles.push_back(glm::uvec3(oddVertexOnEdgeEab, oddVertexOnEdgeEbc, oddVertexOnEdgeEca));
        }

        // after that:
        _triangleIndices = newTriangles;
        _vertexPositions = newVertices;
        recomputePerVertexNormals();
        recomputePerVertexTextureCoordinates();
    }


    void subdivideLoop() {
        // TODO: Implement here the Loop subdivision instead of the straightforward Linear Subdivision.
        // You can have a look at the Linear Subdivision function to take some inspiration from it.
        //
        // A few recommendations / advices (note that the following regards a simple implementation that does not handle boundaries, you can adapt it if you want to handle those):
        // I) start by declaring a vector of new positions "newVertices" and a vector of new triangles "newTriangles".
        //    Do not mix the new quantities and the old ones.
        //    At the end, replace _vertexPositions by newVertices and _triangleIndices by newTriangles, just as it is done in subdivideLinear().
        //    This will help you writing clean code.
        //    Remember: In the Loop subdivision scheme, a new position (in the output mesh at level k+1) is a linear combination of the old vertices positions (at level k).
        //    So, you should NEVER (!!!!!) have in your code something like: newVertices[ v ] += newVertices[ v_neighbor ] * weight;
        // II) Compute the neighbors of all the even vertices. You can use a structure such as "std::vector< std::set< unsigned int > > vertex_neighbors" for example.
        //    This will give you the valence n of a given even vertex v, and the value of the coefficient alpha_n that you need to use in the computation of the new position for v.
        // III) Compute the new positions for the even vertices. If you compute the even vertices first, you will not be tempted to consider the odd vertices as their neighbors (that would be a -- very common, mistake).
        // IV) Process all triangles, insert the odd vertices, compute their position using the subdivision mask, and create four new triangles per old triangle.
        //    You can get inspiration from subdivideLinear() for that part.
        //
        // Good luck! Do not hesitate asking questions, we are here to help you.

        std::vector<glm::vec3> newVertices;
        std::vector<glm::vec3> intermediatVertices;
        std::vector<glm::uvec3> newTriangles;

        std::map< int, std::vector<glm::vec3> > oldVertexOnTriangle;
        for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
            // In this loop we are going to create, to each vertex, an array with the neighboring vertices
            unsigned int a = _triangleIndices[tIt][0];
            unsigned int b = _triangleIndices[tIt][1];
            unsigned int c = _triangleIndices[tIt][2];

            // Create the arrays
            if (oldVertexOnTriangle.find(a) == oldVertexOnTriangle.end()) {
                std::vector<glm::vec3> vertices;
                oldVertexOnTriangle[a] = vertices;
            }

            if (oldVertexOnTriangle.find(b) == oldVertexOnTriangle.end()) {
                std::vector<glm::vec3> vertices;
                oldVertexOnTriangle[b] = vertices;
            }

            if (oldVertexOnTriangle.find(c) == oldVertexOnTriangle.end()) {
                std::vector<glm::vec3> vertices;
                oldVertexOnTriangle[c] = vertices;
            }

            // Check if there's no repeated values 
            if (std::find(oldVertexOnTriangle[a].begin(), oldVertexOnTriangle[a].end(), _vertexPositions[b]) == oldVertexOnTriangle[a].end()) {
                oldVertexOnTriangle[a].push_back(_vertexPositions[b]);
            }
            if (std::find(oldVertexOnTriangle[a].begin(), oldVertexOnTriangle[a].end(), _vertexPositions[c]) == oldVertexOnTriangle[a].end()) {
                oldVertexOnTriangle[a].push_back(_vertexPositions[c]);
            }
            if (std::find(oldVertexOnTriangle[b].begin(), oldVertexOnTriangle[b].end(), _vertexPositions[a]) == oldVertexOnTriangle[b].end()) {
                oldVertexOnTriangle[b].push_back(_vertexPositions[a]);
            }
            if ((std::find(oldVertexOnTriangle[b].begin(), oldVertexOnTriangle[b].end(), _vertexPositions[c])) == oldVertexOnTriangle[b].end()) {
                oldVertexOnTriangle[b].push_back(_vertexPositions[c]);
            }
            if ((std::find(oldVertexOnTriangle[c].begin(), oldVertexOnTriangle[c].end(), _vertexPositions[a])) == oldVertexOnTriangle[c].end()) {
                oldVertexOnTriangle[c].push_back(_vertexPositions[a]);
            }
            if ((std::find(oldVertexOnTriangle[c].begin(), oldVertexOnTriangle[c].end(), _vertexPositions[b])) == oldVertexOnTriangle[c].end()) {
                oldVertexOnTriangle[c].push_back(_vertexPositions[b]);
            }
        }

        // Adjust old vertices -> use the neighboring to calculate the new position
        for (auto const& it : oldVertexOnTriangle) {

            int vertex_index = it.first;
            std::vector<glm::vec3> vertices = it.second;

            float alpha;

            if (vertices.size() > 3) {
                alpha = 3.f / (8.f * (float)vertices.size());
            }
            else {
                alpha = 3.0f / 16.0f;
            }

            glm::vec3 oldV = (1 - ((float)vertices.size()) * alpha) * _vertexPositions[vertex_index];
            glm::vec3 aux = glm::vec3(0, 0, 0);

            for (int i = 0; i < vertices.size(); i++) {
                oldV += alpha * vertices[i];
            }
            newVertices.push_back(oldV);

        }

        struct Edge {
            unsigned int a, b;
            Edge(unsigned int c, unsigned int d) : a(std::min<unsigned int>(c, d)), b(std::max<unsigned int>(c, d)) {}
            bool operator < (Edge const& o) const { return a < o.a || (a == o.a && b < o.b); }
            bool operator == (Edge const& o) const { return a == o.a && b == o.b; }
        };

        std::map< Edge, int > newVertexOnTriangle;
        std::map< Edge, int > oddVertexOnTrianglePoint;
        // Here we are going to create the new odd vertices
        for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
            unsigned int a = _triangleIndices[tIt][0];
            unsigned int b = _triangleIndices[tIt][1];
            unsigned int c = _triangleIndices[tIt][2];

            // Here we store the position of each new vertice
            int oddVertexEab = -1;
            int oddVertexEbc = -1;
            int oddVertexEca = -1;

            // This procedure repats for each edge of the triangle
            Edge Eab(a, b);
            // First we check if we already saw this edge
            if (newVertexOnTriangle.find(Eab) == newVertexOnTriangle.end()) {
                // If we didn't, we are going to add the average to the new vertices (if we don't see it again, it's the boundary case)
                // We map the position of this new Vertex to the Edge so we can acess it later
                // We also map the c (first point of the triangle) to the Edge for the same reason
                newVertices.push_back((_vertexPositions[a] + _vertexPositions[b]) / 2.f);
                oddVertexEab = newVertices.size() - 1;
                newVertexOnTriangle[Eab] = oddVertexEab;
                oddVertexOnTrianglePoint[Eab] = c;
            }
            else {
                // If we already saw this Edge, so we must have already a point in our newVertices array mapped and point of the triangle to calculate
                oddVertexEab = newVertexOnTriangle[Eab]; // old position of oddvertex
                newVertices[oddVertexEab] = (3.f / 8.f) * (_vertexPositions[a] + _vertexPositions[b]) +
                    (1.f / 8.f) * _vertexPositions[oddVertexOnTrianglePoint[Eab]] + (1.f / 8.f) * _vertexPositions[c];
            }


            Edge Ebc(b, c);
            if (newVertexOnTriangle.find(Ebc) == newVertexOnTriangle.end()) {
                newVertices.push_back((_vertexPositions[b] + _vertexPositions[c]) / 2.f);
                oddVertexEbc = newVertices.size() - 1;
                newVertexOnTriangle[Ebc] = oddVertexEbc;
                oddVertexOnTrianglePoint[Ebc] = a;
            }
            else {
                oddVertexEbc = newVertexOnTriangle[Ebc]; // old position of oddvertex
                newVertices[oddVertexEbc] = (3.f / 8.f) * (_vertexPositions[b] + _vertexPositions[c]) +
                    (1.f / 8.f) * _vertexPositions[oddVertexOnTrianglePoint[Ebc]] + (1.f / 8.f) * _vertexPositions[a];
            }

            Edge Eca(c, a);
            if (newVertexOnTriangle.find(Eca) == newVertexOnTriangle.end()) {
                newVertices.push_back((_vertexPositions[c] + _vertexPositions[a]) / 2.f);
                oddVertexEca = newVertices.size() - 1;
                newVertexOnTriangle[Eca] = oddVertexEca;
                oddVertexOnTrianglePoint[Eca] = b;
            }
            else {
                oddVertexEca = newVertexOnTriangle[Eca]; // old position of oddvertex
                newVertices[oddVertexEca] = (3.f / 8.f) * (_vertexPositions[c] + _vertexPositions[a]) +
                    (1.f / 8.f) * _vertexPositions[oddVertexOnTrianglePoint[Eca]] + (1.f / 8.f) * _vertexPositions[b];
            }
            newTriangles.push_back(glm::uvec3(a, oddVertexEab, oddVertexEca));
            newTriangles.push_back(glm::uvec3(oddVertexEab, b, oddVertexEbc));
            newTriangles.push_back(glm::uvec3(oddVertexEca, oddVertexEbc, c));
            newTriangles.push_back(glm::uvec3(oddVertexEab, oddVertexEbc, oddVertexEca));
        }

        _triangleIndices = newTriangles;
        _vertexPositions = newVertices;
        recomputePerVertexNormals();
        recomputePerVertexTextureCoordinates();

    }

private:
    std::vector<glm::vec3> _vertexPositions;
    std::vector<glm::vec3> _vertexNormals;
    std::vector<glm::vec2> _vertexTexCoords;
    std::vector<glm::uvec3> _triangleIndices;

    GLuint _vao = 0;
    GLuint _posVbo = 0;
    GLuint _normalVbo = 0;
    GLuint _texCoordVbo = 0;
    GLuint _ibo = 0;
};

// utility: loader
void loadOFF(const std::string& filename, std::shared_ptr<Mesh> meshPtr);

#endif  // MESH_H
