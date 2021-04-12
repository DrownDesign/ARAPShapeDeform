#include <cstdio>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <string>
#include <glm/glm.hpp>
#include <iostream>


#include "Mesh.h"

using namespace std;
using namespace glm;

//Loaded Mesh Data
vector<vec3> loadedVertices, loadedNormals;
vector<vec3> outVertices, outNorms;
vector<int> loadedVertexIndices, loadedUVInidces, loadedNormalIndices;
vector<vec2> tempUVs;



bool loadModel(const char * path) {

	printf("Importing Model...\n");	

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
			loadedVertices.push_back(vertex);
		}
		/*else if (strcmp(lineHeader, "vt") == 0) {
			vec2 uv;
			fscanf(file, "%f %f\n", &uv.x, &uv.y);
			tempUVs.push_back(uv);
		}*/
		else if (strcmp(lineHeader, "vn") == 0) {
			vec3 norm;
			fscanf(file, "%f %f %f\n", &norm.x, &norm.y, &norm.z);
			loadedNormals.push_back(norm);
		}
		else if (strcmp(lineHeader, "f") == 0) {
			string vertex1, vertex2, vertex3;
			unsigned int vertexIndex[3], uvIndex[3], normalIndex[3];
			int matches = fscanf(file, "%d//%d %d//%d %d//%d\n", &vertexIndex[0], &normalIndex[0],
				&vertexIndex[1], &normalIndex[1],
				&vertexIndex[2], &normalIndex[2]);
			
			loadedVertexIndices.push_back(vertexIndex[0]);
			loadedVertexIndices.push_back(vertexIndex[1]);
			loadedVertexIndices.push_back(vertexIndex[2]);

			loadedNormalIndices.push_back(normalIndex[0]);
			loadedNormalIndices.push_back(normalIndex[1]);
			loadedNormalIndices.push_back(normalIndex[2]);
		}

		for (unsigned int i = 0; i < loadedVertexIndices.size(); i++) {
			int vertexIndex = loadedVertexIndices[i];
			vec3 vertex = loadedVertices[vertexIndex - 1];
			outVertices.push_back(vertex);
		}

		for (unsigned int i = 0; i < loadedNormalIndices.size(); i++) {
			int normalIndex = loadedNormalIndices[i];
			vec3 normal = loadedNormals[normalIndex - 1];
			outNorms.push_back(normal);
		}


	}
	fclose(file);
	printf("End of Load\n");
	return true;
}

//Core code for assignment
vector<vec3> deformedVertices, deformedNormals;

void arapDeform() {
	
	//Testing export logic 
	deformedVertices = loadedVertices;
	deformedNormals = loadedNormals;

}

//Exports the model after deformation
bool exportModel() {
	printf("Exporting...\n");

	string path = "./ExportObj/DeformedMesh.obj";
	
	ofstream objExport;
	objExport.open(path);

	//Header
	objExport << "####\n";
	objExport << "#\n";
	objExport << "#Resulting Export of ARAP Deformation\n";
	objExport << "#\n";
	objExport << "####\n";
	objExport << "# Object DeformedMesh.obj\n";
	objExport << "#\n";
	objExport << "# Vertices: %d\n";
	objExport << "#\n";
	objExport << "# Faces: %d\n";
	objExport << "#\n";
	objExport << "####\n";

	for (int i = 0; i < loadedVertices.size(); i++) {
		vec3 v = loadedVertices.at(i);
		objExport << "v " + to_string(v.x) + " " + to_string(v.y) + " " + to_string(v.z) + "\n";
	}

	for (int i = 0; i < loadedNormals.size(); i++) {
		vec3 v = loadedNormals.at(i);
		objExport << "vn " + to_string(v.x) + " " + to_string(v.y) + " " + to_string(v.z) + "\n";
	}

	for (int i = 0; i < loadedVertexIndices.size(); i+=3) {
		vec3 vertIndex(loadedVertexIndices.at(i), loadedVertexIndices.at(i+1), loadedVertexIndices.at(i+2));
		vec3 normIndex(loadedNormalIndices.at(i), loadedNormalIndices.at(i+1), loadedNormalIndices.at(i+2));
	
		objExport << "f " + to_string(vertIndex.x) + "//" + to_string(normIndex.x) + " " + 
			to_string(vertIndex.y) + "//" + to_string(normIndex.y) + " " + 
			to_string(vertIndex.z) + "//" + to_string(normIndex.z) + "\n";
	}

	//Start file here

	objExport.close();


	printf("Successful Export");
	return true;
}

//Runs on startup
int main(int argc, char **argv) {

	const char *filename = "./InputObj/bar1.obj";

	//Load Model
	loadModel(filename);

	//Deform Model
	arapDeform();

	//Export Model
	exportModel();

	exit(0);
}