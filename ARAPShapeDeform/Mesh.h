#pragma once
#include <vector>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

class Mesh
{
	
public:
	vector<vec3> vertices;
	vector<vec3> indices;
	vector<vec2> texCoords;
	vector<vec3> faces;

	Mesh(vector<vec3> vertices, vector<vec3> indices, vector<vec2> texCoords, vector<vec3> faces);
};

