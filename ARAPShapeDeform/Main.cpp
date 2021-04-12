#include <cstdio>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <string>
#include <glm/glm.hpp>


#include "Mesh.h"

using namespace std;
using namespace glm;

//Mesh m_mesh;

bool loadModel(const char * path) {

	printf("Importing Model...\n");

	vector<int> vertexIndices, uvIndices, normalIndices;
	vector<vec3> tempVertices, tempNormals;
	vector<vec2> tempUVs;

	vector<vec3> outVertices, outNorms;

	FILE * file = fopen(path, "r");
	if (file == NULL) {
		printf("Error Opening File\n");
		return false;
	}

	while (1) {
		char lineHeader[128];
		int res = fscanf(file, "%s", lineHeader);
		
		if (res == EOF) {
			break;
		}

		if (strcmp(lineHeader, "v") == 0) {
			vec3 vertex;
			fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
			tempVertices.push_back(vertex);
		}
		else if (strcmp(lineHeader, "vt") == 0) {
			vec2 uv;
			fscanf(file, "%f %f\n", &uv.x, &uv.y);
			tempUVs.push_back(uv);
		}
		else if (strcmp(lineHeader, "vn") == 0) {
			vec3 norm;
			fscanf(file, "%f %f %f\n", &norm.x, &norm.y, &norm.z);
			tempNormals.push_back(norm);
		}
		else if (strcmp(lineHeader, "f") == 0) {
			string vertex1, vertex2, vertex3;
			unsigned int vertexIndex[3], uvIndex[3], normalIndex[3];
			int matches = fscanf(file, "%d//%d %d//%d %d//%d\n", &vertexIndex[0], &normalIndex[0],
				&vertexIndex[1], &normalIndex[1],
				&vertexIndex[2], &normalIndex[2]);
			
			vertexIndices.push_back(vertexIndex[0]);
			vertexIndices.push_back(vertexIndex[1]);
			vertexIndices.push_back(vertexIndex[2]);

			normalIndices.push_back(normalIndex[0]);
			normalIndices.push_back(normalIndex[1]);
			normalIndices.push_back(normalIndex[2]);
		}

		for (unsigned int i = 0; i < vertexIndices.size(); i++) {
			int vertexIndex = vertexIndices[i];
			vec3 vertex = tempVertices[vertexIndex - 1];
			outVertices.push_back(vertex);
		}

		for (unsigned int i = 0; i < normalIndices.size(); i++) {
			int normalIndex = normalIndices[i];
			vec3 normal = tempNormals[normalIndex - 1];
			outNorms.push_back(normal);
		}

		//m_mesh = Mesh(outVertices, outNorms);

	}
	fclose(file);
	printf("End of Load\n");
	return true;
}

//Exports the model after deformation
bool exportModel() {
	printf("Exporting...\n");

	const char *path = "./ExportObj/DeformedMesh.obj";

	
	return true;
}

//Runs on startup
int main(int argc, char **argv) {

	const char *filename = "./InputObj/bar1.obj";

	//Load Model
	loadModel(filename);

	//Deform Model

	//Export Model
	exportModel();

	exit(0);
}