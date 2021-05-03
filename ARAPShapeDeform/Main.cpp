#include <cstdio>
#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <glm/glm.hpp>
#include <iostream>

//Eigen
#include <Eigen/Dense>

#include "Mesh.h"

using namespace std;
using namespace glm;
using namespace Eigen;

//Loaded Mesh Data
vector<vec3> loadedVertices, loadedNormals;
vector<vec3> outVertices, outNorms;
vector<int> loadedVertexIndices, loadedUVInidces, loadedNormalIndices;
map<int, set<int>> edges;
vector<vec2> tempUVs;

bool loadModel(const char * path) {

	printf("Importing Model...\n");

	FILE * file = fopen(path, "r");
	if (file == NULL) {
		printf("Error Opening Mesh\n");
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
			int vertexIndex[3], normalIndex[3];
			int matches = fscanf(file, "%d//%d %d//%d %d//%d\n", &vertexIndex[0], &normalIndex[0],
				&vertexIndex[1], &normalIndex[1],
				&vertexIndex[2], &normalIndex[2]);

			loadedVertexIndices.push_back(vertexIndex[0]);
			loadedVertexIndices.push_back(vertexIndex[1]);
			loadedVertexIndices.push_back(vertexIndex[2]);

			loadedNormalIndices.push_back(normalIndex[0]);
			loadedNormalIndices.push_back(normalIndex[1]);
			loadedNormalIndices.push_back(normalIndex[2]);

			//Make map of adj vertices for adj matrix
			edges[vertexIndex[0]-1].insert(vertexIndex[1]-1);
			edges[vertexIndex[0]-1].insert(vertexIndex[2]-1);

			edges[vertexIndex[1]-1].insert(vertexIndex[0]-1);
			edges[vertexIndex[1]-1].insert(vertexIndex[2]-1);

			edges[vertexIndex[2]-1].insert(vertexIndex[0]-1);
			edges[vertexIndex[2]-1].insert(vertexIndex[1]-1);
		}
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

	fclose(file);
	printf("End of Load\n");
	return true;
}

//Loaded Constraints vector
vector<int> vertTypes;

bool loadVertexTypes(const char *filepath) {

	printf("Loading Constraints...\n");

	//Check file
	FILE * file = fopen(filepath, "r");
	if (file == NULL) {
		printf("Error Opening Constraints\n");
		return false;
	}

	while (1) {
		char lineHeader[128];

		int res = fscanf(file, "%s", lineHeader);

		if (res == EOF) {
			break;
		}

		if (strcmp(lineHeader, "#") == 0) {
			//Skip comments
			char buffer[100];
			fgets(buffer, 100, file);
			//cout << buffer << endl;
		}
		else {
			int value, v1;
			v1 = stof(lineHeader);
			fscanf(file, "%i", &value);

			//Hacky fix for line skipping but it works
			vertTypes.push_back(v1);
			vertTypes.push_back(value);
		}

	}

	return true;
}

//Loaded Deformation Matrix
MatrixXf deformationMatrix(4, 4);

bool loadDeformation(const char *filepath) {

	printf("Loading Deformation...\n");

	//Check file
	FILE * file = fopen(filepath, "r");
	if (file == NULL) {
		printf("Error Opening Constraints\n");
		return false;
	}

	int count = 0;

	while (1) {
		char lineHeader[128];

		int res = fscanf(file, "%s", lineHeader);

		if (res == EOF) {
			break;
		}

		if (strcmp(lineHeader, "#") == 0) {
			//Skip comments
			char buffer[100];
			fgets(buffer, 100, file);
			//cout << buffer << endl;
		}
		else {
			float v1, v2, v3, v4;
			v1 = stof(lineHeader);
			fscanf(file, "%f %f %f\n", &v2, &v3, &v4);

			//Load values to matrix
			deformationMatrix(count, 0) = v1;
			deformationMatrix(count, 1) = v2;
			deformationMatrix(count, 2) = v3;
			deformationMatrix(count, 3) = v4;

			//printf("%f %f %f %f\n", v1, v2, v3, v4);

			count++;

		}
	}

	fclose(file);

	//Print Matrix
	//cout << deformationMatrix << endl;
	//printf("Matrix Loaded\n");


	return true;
}

//Laplacian Matrix
MatrixXf laplacianMatrix;

//Generates a Laplacian matrix for the loaded mesh
void generateLaplacian() {
	
	int numVertices = loadedVertices.size();

	//Adjacency Matrix <- seems to work?
	MatrixXi adjMatrix(numVertices, numVertices);
	adjMatrix.fill(0);

	for (int j = 0; j < numVertices ; j++) {
		for (int i = 0; i < numVertices ; i++) {
			bool edgeExists = false;

			auto neighbours = edges.find(i);

			if (neighbours != edges.end()) {
				set<int> test = neighbours->second;
				for (int smt : test) {
					if (smt == j) {
						edgeExists = true;
						continue;
					}
				}
			}

			/*for (const auto& neighbours : edges) {
				if (neighbours.first == i) {
					set<int> test = neighbours.second;
					for (int smt : test) {
						if (smt == j) {
							edgeExists = true;
							continue;
						}
					}
				}
			}*/

			if (edgeExists) {
				adjMatrix(i, j) = 1;
			}
			else {
				adjMatrix(i, j) = 0;
			}

		}
	}

	printf("\nAdjacency Matrix:\n");
	cout << adjMatrix << endl;
	printf("-------------------------------------------\n");


	//Degree Matrix
	MatrixXi degMatrix(numVertices, numVertices);
	degMatrix.fill(0);
	
	for (int j = 0; j < numVertices; j++) {
		for (int i = 0; i < numVertices; i++) {
			
			int degree = 0;

			if (i==j) {

				for (const auto &nbrs : edges) {
					if (nbrs.first == i) {
						degree = nbrs.second.size();
					}
				}
				
				degMatrix(i, j) = degree;
			}

		}
	}

	printf("\nDegree Matrix (from adjacency):\n");
	cout << degMatrix << endl;
	printf("-------------------------------------------\n");

	laplacianMatrix = degMatrix.cast<float>() - adjMatrix.cast<float>();

	printf("\nLaplacian Matrix (L = D - A):\n");
	cout << laplacianMatrix << endl;
	printf("-------------------------------------------\n");

}

vector<vec3> diffCoords;
VectorXf xDiffPositions, yDiffPositions, zDiffPositions;

void generateDifCoords() {

	//Seperate positions into x,y,z vectors
	VectorXf xPositions(loadedVertices.size()), yPositions(loadedVertices.size()), zPositions(loadedVertices.size());
	
	for (int i = 0; i < loadedVertices.size(); i++) {
		vec3 originalPosition = loadedVertices.at(i);

		xPositions(i) = originalPosition.x;
		yPositions(i) = originalPosition.y;
		zPositions(i) = originalPosition.z;

	}

	//Calculate differentials
	xDiffPositions = laplacianMatrix * xPositions;
	yDiffPositions = laplacianMatrix * yPositions;
	zDiffPositions = laplacianMatrix * zPositions;

	for (int i = 0; i < loadedVertices.size(); i++) {
		vec3 diffPosition(xDiffPositions(i), yDiffPositions(i), zDiffPositions(i));
		
		diffCoords.push_back(diffPosition);

		
		/*printf("Original Position: %f %f %f		New Position: %f %f %f\n", loadedVertices.at(i).x, loadedVertices.at(i).y, loadedVertices.at(i).z, 
			diffCoords.at(i).x, diffCoords.at(i).y, diffCoords.at(i).z);*/
	}
	
}

//Find optimal rotation (local step)
void findRotation() {

}

vector<vec3> recoveredPositions;

//Recover the surface of the model by solving a linear system
void recoverSurface() {

	//Find Constraints
	vector<int> constraintIndices;
	int index = 0;
	for (int v : vertTypes) {
		if (v == 0) { constraintIndices.push_back(index); }
		index++;
	}

	MatrixXf tempLaplacian = laplacianMatrix;

	//Alter Laplacian to match
	for (int i : constraintIndices) {
		for (int x = 0; x < loadedVertices.size(); x++) {
			tempLaplacian(i, x) = 0;

			if (x == i) {
				tempLaplacian(i, x) = 1;
			}
		}
	}

	//Invert Laplacian
	MatrixXf InvLaplacian = tempLaplacian.inverse();

	printf("\nTemp Laplacian:\n");
	cout << tempLaplacian << endl;
	printf("-------------------------------------------\n");

	printf("\nInverse Laplacian:\n");
	cout << InvLaplacian << endl;
	printf("-------------------------------------------\n");

	//Fixed Example
	xDiffPositions(2) = 15;
	xDiffPositions(5) = 20;
	yDiffPositions(2) = 15;
	yDiffPositions(5) = 30;
	zDiffPositions(2) = 3;
	zDiffPositions(5) = 3;


	cout << xDiffPositions << endl;

	//Multiply LT by d(x,y,z)
	VectorXf recoveredX = InvLaplacian * xDiffPositions;
	VectorXf recoveredY = InvLaplacian * yDiffPositions;
	VectorXf recoveredZ = InvLaplacian * zDiffPositions;

	cout << recoveredX << endl;

	for (int i = 0; i < loadedVertices.size(); i++) {
		vec3 recoveredPos(recoveredX(i), recoveredY(i), recoveredZ(i));
		recoveredPositions.push_back(recoveredPos);
		printf("Index: %i, Position: %f %f %f\n", i, recoveredPos.x, recoveredPos.y, recoveredPos.z);
	}
}

//Exports the model after deformation
bool exportModel() {
	printf("Exporting...\n");

	string path = "./ExportObj/4Face.obj";

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

	for (int i = 0; i < recoveredPositions.size(); i++) {
		vec3 v = recoveredPositions.at(i);
		objExport << "v " + to_string(v.x) + " " + to_string(v.y) + " " + to_string(v.z) + "\n";
	}

	for (int i = 0; i < loadedNormals.size(); i++) {
		vec3 v = loadedNormals.at(i);
		objExport << "vn " + to_string(v.x) + " " + to_string(v.y) + " " + to_string(v.z) + "\n";
	}

	for (int i = 0; i < loadedVertexIndices.size(); i += 3) {
		vec3 vertIndex(loadedVertexIndices.at(i), loadedVertexIndices.at(i + 1), loadedVertexIndices.at(i + 2));
		vec3 normIndex(loadedNormalIndices.at(i), loadedNormalIndices.at(i + 1), loadedNormalIndices.at(i + 2));

		objExport << "f " + to_string(vertIndex.x) + "//" + to_string(normIndex.x) + " " +
			to_string(vertIndex.y) + "//" + to_string(normIndex.y) + " " +
			to_string(vertIndex.z) + "//" + to_string(normIndex.z) + "\n";
	}

	//Start file here

	objExport.close();


	printf("Successful Export\n");
	return true;
}

//Runs on startup
int main(int argc, char **argv) {

	const char *meshFile = "./InputObj/4Face.obj";
	const char *vTypesFile = "./InputObj/4Face.sel";
	const char *deformationFile = "./InputObj/4Face.def";

	//Load Files
	loadModel(meshFile);
	loadVertexTypes(vTypesFile);
	loadDeformation(deformationFile);

	//Generate Laplacian Matrix
	generateLaplacian();

	//Generate Differential Coordinates
	generateDifCoords();

	//Find optimal rotation (Local) Need to put in loop with global
	findRotation();

	//Recover surface (Global)
	recoverSurface();

	//Export Model
	exportModel();

	exit(0);
}