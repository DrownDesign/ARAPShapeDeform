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

using namespace std;
using namespace glm;
using namespace Eigen;

//SETTINGS
const bool GLOBAL_STEP = true;
const int OFFSET = 3;
const char *meshFile = "./InputObj/4Face.obj";
const char *vTypesFile = "./InputObj/4Face.sel";
const char *deformationFile = "./InputObj/4Face.def";
const int iterations = 1;
const int DEBUG_ROT = true;



//Loaded Mesh Data
vector<vec3> loadedVertices, loadedVerticesPrime, loadedNormals;
vector<int> loadedVertexIndices, loadedNormalIndices;
map<int, set<int>> neighbours;
vector<MatrixXf> rotations;
Matrix4f initialRotation;

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
			neighbours[vertexIndex[0]-1].insert(vertexIndex[1]-1);
			neighbours[vertexIndex[0]-1].insert(vertexIndex[2]-1);

			neighbours[vertexIndex[1]-1].insert(vertexIndex[0]-1);
			neighbours[vertexIndex[1]-1].insert(vertexIndex[2]-1);

			neighbours[vertexIndex[2]-1].insert(vertexIndex[0]-1);
			neighbours[vertexIndex[2]-1].insert(vertexIndex[1]-1);
		}
	}

	fclose(file);

	//Same values at start
	loadedVerticesPrime = loadedVertices;
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
			count++;

		}
	}

	fclose(file);
	
	//For debugging
	initialRotation <<
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;

	for (int i = 0; i < loadedVertices.size(); i++) {
		if (DEBUG_ROT) {
			rotations.push_back(initialRotation);

		}
		else {
			rotations.push_back(deformationMatrix);
		}
	}
	
	/*printf("\nDeformation Matrix:\n");
	cout << deformationMatrix << endl;
	printf("-------------------------------------------\n");*/
	return true;
}

MatrixXf bMatrix;

void intializeBMatrix() {
	printf("Initializing B Matrix\n");
	bMatrix.resize(loadedVertices.size(), 3);
	bMatrix.fill(0);

	for (int i = 0; i < loadedVertices.size(); i++) {
		auto pair = neighbours.find(i);

		if (pair != neighbours.end()) {
			set<int> nSet = pair->second;
			for (int vNIndex : nSet) {
				Vector4f current(loadedVertices.at(i).x, loadedVertices.at(i).y, loadedVertices.at(i).z,1);
				Vector4f vN(loadedVertices.at(vNIndex).x, loadedVertices.at(vNIndex).y, loadedVertices.at(vNIndex).z,1);
				Vector4f newVec(0.5f * (rotations.at(i) + rotations.at(vNIndex)) * (current - vN));
				//Check that this should be  += not =
				bMatrix(i, 0) += newVec(0, 0);
				bMatrix(i, 1) += newVec(1, 0);
				bMatrix(i, 2) += newVec(2, 0);
			}
		}
	}

	/*printf("\nB Matrix:\n");
	cout << bMatrix << endl;
	printf("-------------------------------------------\n");*/
}

MatrixXf laplacianMatrix;

//Generates a Laplacian matrix for the loaded mesh
void generateLaplacian() {
	printf("Building Laplacian\n");
	int numVertices = loadedVertices.size();

	//Adjacency Matrix <- seems to work?
	MatrixXi adjMatrix(numVertices, numVertices);
	adjMatrix.fill(0);

	for (int j = 0; j < numVertices ; j++) {
		for (int i = 0; i < numVertices ; i++) {
			bool edgeExists = false;

			auto vN = neighbours.find(i);

			if (vN != neighbours.end()) {
				set<int> test = vN->second;
				for (int smt : test) {
					if (smt == j) {
						edgeExists = true;
						continue;
					}
				}
			}

			if (edgeExists) {
				adjMatrix(i, j) = 1;
			}
			else {
				adjMatrix(i, j) = 0;
			}

		}
	}

	/*printf("\nAdjacency Matrix:\n");
	cout << adjMatrix << endl;
	printf("-------------------------------------------\n");*/


	//Degree Matrix
	MatrixXi degMatrix(numVertices, numVertices);
	degMatrix.fill(0);
	
	for (int j = 0; j < numVertices; j++) {
		for (int i = 0; i < numVertices; i++) {
			
			int degree = 0;

			if (i==j) {

				for (const auto &nbrs : neighbours) {
					if (nbrs.first == i) {
						degree = nbrs.second.size();
					}
				}
				
				degMatrix(i, j) = degree;
			}

		}
	}

	/*printf("\nDegree Matrix (from adjacency):\n");
	cout << degMatrix << endl;
	printf("-------------------------------------------\n");*/

	laplacianMatrix = degMatrix.cast<float>() - adjMatrix.cast<float>();

	/*printf("\nLaplacian Matrix (L = D - A):\n");
	cout << laplacianMatrix << endl;
	printf("-------------------------------------------\n");*/

}

MatrixXf Lc, Bc;

void generateGlobalConstraints() {
	printf("Generating Constraints\n");
	//Find Constraints
	vector<int> constraintIndices;
	int index = 0;
	for (int v : vertTypes) {
		if (v == 0) { constraintIndices.push_back(index); }
		index++;
	}

	//Constrain v wrt laplacian
	Lc = laplacianMatrix;
	for (int i : constraintIndices) {
		for (int x = 0; x < loadedVertices.size(); x++) {
			Lc(i, x) = 0;

			if (x == i) {
				Lc(i, x) = 1;
			}
		}
	}

	/*printf("\nLc Matrix:\n");
	cout << Lc << endl;
	printf("-------------------------------------------\n");*/


	//Constrain v wrt BMatrix
	Bc = bMatrix;

	for (int i : constraintIndices) {
		for (int x = 0; x < Bc.cols(); x++) {
			Bc(i, x) = OFFSET + i;
		}
	}
}

void generateConstraints() {
	printf("Generating Constraints\n");
	//Find Constraints
	vector<int> constraintIndices;
	int index = 0;
	for (int v : vertTypes) {
		if (v == 0) { constraintIndices.push_back(index); }
		index++;
	}

	//Constrain v wrt laplacian
	Lc = laplacianMatrix;
	for (int i : constraintIndices) {
		for (int x = 0; x < loadedVertices.size(); x++) {
			Lc(i, x) = 0;

			if (x == i) {
				Lc(i, x) = 1;
			}
		}
	}

	/*printf("\nLc Matrix:\n");
	cout << Lc << endl;
	printf("-------------------------------------------\n");*/

	//vector<Vector4f> cRotated;
	//for (int i : constraintIndices) {
	//	Vector4f rotatedConstraint(loadedVertices.at(i).x, loadedVertices.at(i).y, loadedVertices.at(i).z, 1);
	//	rotatedConstraint = initialRotation * rotatedConstraint;
	//	cRotated.push_back(rotatedConstraint);
	//}

	////Constrain v wrt BMatrix
	//Bc = bMatrix;

	//for (int i : constraintIndices) {
	//	Vector4f value = cRotated.at(0);
	//	Bc(i, 0) = value(0);
	//	Bc(i, 1) = value(1);
	//	Bc(i, 2) = value(2);
	//}

	//Constrain v wrt BMatrix
	Bc = bMatrix;

	for (int i : constraintIndices) {			
		Bc(i,0) = loadedVertices.at(i).x;
		Bc(i,1) = loadedVertices.at(i).y;
		Bc(i,2) = loadedVertices.at(i).z;
	}

	/*printf("\nBc Matrix:\n");
	cout << Bc << endl;
	printf("-------------------------------------------\n");*/
}

MatrixXf Dc;
vector<Vector3f> vPrimes;

void solvePositions() {
	printf("Solving Positions\n");
	Dc = Lc.inverse() * Bc;

	/*printf("\nDc Matrix:\n");
	cout << Dc << endl;
	printf("-------------------------------------------\n");*/

	for (int i = 0; i < Dc.rows(); i++) {
		Vector3f newVec;
		for (int j = 0; j < Dc.cols(); j++) {
			newVec(j) = Dc(i, j);
		}
		vPrimes.push_back(newVec);
	}
}

void solveRotations() {
	printf("Solving Rotations\n");

	for (int i = 0; i < loadedVertices.size(); i++) {
		vec3 oldV = loadedVertices.at(i);
		Vector3f v(loadedVertices.at(i).x, loadedVertices.at(i).y, loadedVertices.at(i).z);
		Vector3f v1 = vPrimes.at(i);

		//Calculate centroids
		Vector3f vHat, primeHat;
		vHat.fill(0);
		primeHat.fill(0);

		auto pair = neighbours.find(i);
		if (pair != neighbours.end()) {
			set<int> nSet = pair->second;
			for (int vNIndex : nSet) {
				
				//Get v values
				Vector3f oldN(loadedVertices.at(vNIndex).x, loadedVertices.at(vNIndex).y, loadedVertices.at(vNIndex).z);
				vHat += oldN;

				//Get vPrime values
				primeHat += vPrimes.at(vNIndex);

			}
			vHat /= pair->second.size();
			primeHat /= pair->second.size();
		}

		Vector3f cV = v - vHat;
		Vector3f cVPrime = v1 - primeHat;

		MatrixXf covariance = cV * cVPrime.transpose();

		JacobiSVD<MatrixXf> svd(covariance, ComputeThinU | ComputeThinV);
		MatrixXf uMatrix = svd.matrixU();
		MatrixXf vMatrix = svd.matrixV();

		MatrixXf optimalRot = vMatrix * uMatrix.transpose();
		rotations.at(i) = optimalRot;
		
		/*printf("\nOptimal Rotation Matrix:\n");
		cout << optimalRot << endl;
		printf("-------------------------------------------\n");*/

	}

}

//Exports the model after deformation
bool exportModel() {
	printf("Exporting...\n");

	string path = "./ExportObj/deformedMesh.obj";

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

	for (int i = 0; i < Dc.rows(); i++) {
		
		float vX, vY, vZ;

		for (int j = 0; j < Dc.cols(); j++) {
			if (j == 0) {
				vX = Dc(i, 0);
			}
			else if (j == 1) {
				vY = Dc(i, 1);
			}
			else {
				vZ = Dc(i, 2);
			}
		}

		objExport << "v " + to_string(vX) + " " + to_string(vY) + " " + to_string(vZ) + "\n";
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

	//Load Files
	loadModel(meshFile);
	loadVertexTypes(vTypesFile);
	loadDeformation(deformationFile);

	intializeBMatrix();
	generateLaplacian();

	if (GLOBAL_STEP) {
		//Constrain vertices
		generateGlobalConstraints();

		//solve for best positions
		solvePositions();
	}
	else {
		for (int i = 0; i < iterations; i++) {
			printf("Iteration: %i\n", i);
			//Constrain vertices
			generateConstraints();

			//solve for best positions
			solvePositions();

			//solve for best rotations
			solveRotations();
		}
	}
	
	//Export Model
	exportModel();

	exit(0);
}