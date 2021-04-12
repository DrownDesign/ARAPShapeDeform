#include "Mesh.h"

Mesh::Mesh(vector<vec3> vertices, vector<vec3> indices, vector<vec2> texCoords, vector<vec3> faces)
{
	this->vertices = vertices;
	this->indices = indices;
	this->texCoords = texCoords;
	this->faces = faces;
}
