#include <cstdio>
#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <iostream>

//Eigen
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//SETTINGS
const bool GLOBAL_STEP = false;
const int OFFSET = 3;
const char *meshFile = "./InputObj/tri.obj";
const char *vTypesFile = "./InputObj/tri.sel";
const char *deformationFile = "./InputObj/tri.def";
const int iterations = 10;
const int DEBUG_ROT = false;

//Loaded Mesh Data
vector<Vector3f> loadedVertices, loadedNormals;
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
			float x, y, z;
			fscanf(file, "%f %f %f\n", &x, &y, &z);
			loadedVertices.push_back(Vector3f(x, y, z));
		}
		else if (strcmp(lineHeader, "vn") == 0) {
			float xn, yn, zn;
			fscanf(file, "%f %f %f\n", &xn, &yn, &zn);
			loadedNormals.push_back(Vector3f(xn, yn, zn));
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
			neighbours[vertexIndex[0] - 1].insert(vertexIndex[1] - 1);
			neighbours[vertexIndex[0] - 1].insert(vertexIndex[2] - 1);

			neighbours[vertexIndex[1] - 1].insert(vertexIndex[0] - 1);
			neighbours[vertexIndex[1] - 1].insert(vertexIndex[2] - 1);

			neighbours[vertexIndex[2] - 1].insert(vertexIndex[0] - 1);
			neighbours[vertexIndex[2] - 1].insert(vertexIndex[1] - 1);
		}
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

	char line[1024];

	while (fgets(line, sizeof line, file)) {
		int num;
		if (sscanf(line, "%i", &num) == 1) {
			
			vertTypes.push_back(num);
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
				Vector3f current = loadedVertices.at(i);
				Vector4f v4(current(0), current(1), current(2), 1);

				Vector3f neighbour = loadedVertices.at(vNIndex);
				Vector4f n4(neighbour(0), neighbour(1), neighbour(2), 1);


				Vector4f newVec(0.5f * (rotations.at(i) + rotations.at(vNIndex)) * (v4 - n4));
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

	for (int j = 0; j < numVertices; j++) {
		for (int i = 0; i < numVertices; i++) {
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

			if (i == j) {

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
		for (int x = 0; x < Lc.cols(); x++) {
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
		Vector3f current = loadedVertices.at(i);
		Bc(i, 0) = current(0);
		Bc(i, 1) = current(1);
		Bc(i, 2) = current(2);
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

	/*Vector3f a1(0, 0, 0);
	Vector3f a2(1, 0, 0);
	Vector3f a3(0, 1, 0);
	Vector3f b1(0, 0, 0);
	Vector3f b2(0, 1, 0);
	Vector3f b3(-1, 0, 0);
	Matrix3f vHat;
	
	vHat.col(0) = a1; vHat.col(1) = a2; vHat.col(2) = a3;
	Matrix3f vPrime; vPrime.col(0) = b1; vPrime.col(1) = b2; vPrime.col(2) = b3;*/

	for (int i = 0; i < loadedVertices.size(); i++) {

		auto search = neighbours.find(i);
		if (search != neighbours.end()) {
			set<int> nForI = search->second;
			int sizeOfN = nForI.size();

			MatrixXf vHat(3, sizeOfN), primeHat(3, sizeOfN);

			Vector3f currentV = loadedVertices.at(i);
			int nn = 0;
			for (int index : nForI) {
				Vector3f neighbourVertex = loadedVertices.at(index);
				Vector3f e = currentV - neighbourVertex;
				vHat.col(nn) = e;
				nn++;
			}

			nn = 0;
			Vector3f currentPrimeV = vPrimes.at(i);
			for (int index : nForI) {
				Vector3f neighbourVertex = vPrimes.at(index);
				Vector3f e = currentPrimeV - neighbourVertex;
				primeHat.col(nn) = e;
				nn++;
			}
			/*Vector3f a1(0, 0, 0);
			Vector3f a2(1, 0, 0);
			Vector3f a3(0, 1, 0);
			Vector3f b1(0, 0, 0);
			Vector3f b2(0, 1, 0);
			Vector3f b3(-1, 0, 0);
			Matrix3f vHat; vHat.col(0) = a1; vHat.col(1) = a2; vHat.col(2) = a3;
			Matrix3f vPrime; vPrime.col(0) = b1; vPrime.col(1) = b2; vPrime.col(2) = b3;*/

			cout << "vHat: " << vHat << endl;
			cout << "pHat: " << primeHat << endl;


			MatrixXf covariance = vHat * primeHat.transpose();
			JacobiSVD<MatrixXf> svd(covariance, ComputeThinU | ComputeThinV);
			MatrixXf uMatrix = svd.matrixU(), vMatrix = svd.matrixV();
			MatrixXf optimalRot = vMatrix * uMatrix.transpose();
			cout << optimalRot << endl;
			rotations.at(i) = optimalRot;
		}

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

	for (int i = 0; i < vPrimes.size(); i++) {
		Vector3f v = vPrimes.at(i);
		objExport << "v " + to_string(v(0)) + " " + to_string(v(1)) + " " + to_string(v(2)) + "\n";
	}

	for (int i = 0; i < loadedNormals.size(); i++) {
		Vector3f v = loadedNormals.at(i);
		objExport << "vn " + to_string(v(0)) + " " + to_string(v(1)) + " " + to_string(v(2)) + "\n";
	}

	for (int i = 0; i < loadedVertexIndices.size(); i += 3) {
		objExport << "f " + to_string(loadedVertexIndices.at(i)) + "//" + to_string(loadedNormalIndices.at(i)) + " " +
			to_string(loadedVertexIndices.at(i + 1)) + "//" + to_string(loadedNormalIndices.at(i + 1)) + " " +
			to_string(loadedVertexIndices.at(i + 2)) + "//" + to_string(loadedNormalIndices.at(i + 2)) + "\n";
	}

	objExport.close();

	printf("Successful Export\n");
	return true;
}

void runGlobalTest() {
	//Load Files
	loadModel(meshFile);
	loadVertexTypes(vTypesFile);
	loadDeformation(deformationFile);

	intializeBMatrix();
	generateLaplacian();

	//Constrain vertices
	generateGlobalConstraints();

	//solve for best positions
	solvePositions();

	//Export Model
	exportModel();
}

void runARAP() {
	//Load Files
	loadModel(meshFile);
	loadVertexTypes(vTypesFile);
	loadDeformation(deformationFile);

	intializeBMatrix();
	generateLaplacian();

	for (int i = 0; i < iterations; i++) {
		printf("Iteration: %i\n", i);
		//Constrain vertices
		generateConstraints();

		//solve for best positions
		solvePositions();

		//solve for best rotations
		solveRotations();
		loadedVertices = vPrimes;
	}

	//Export Model
	exportModel();
}

//Runs on startup
int main(int argc, char **argv) {

	if (GLOBAL_STEP) {
		runGlobalTest();
	}
	else {
		runARAP();
	}

	exit(0);
}