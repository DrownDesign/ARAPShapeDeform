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
const int OFFSET = 3;
const char *meshFile = "./InputObj/cylinder.obj";
const char *vTypesFile = "./InputObj/cylinder.sel";
const char *deformationFile = "./InputObj/cylinder.def";
const int iterations = 10;
const int DEBUG_ROT = false;

//Loaded Mesh Data
vector<Vector3f> loadedVertices, loadedNormals;
vector<int> loadedVertexIndices, loadedNormalIndices;
map<int, set<int>> neighbours;
vector<MatrixXf> rotations;
Matrix3f initialRotation;
vector<Vector3f> vPrimes;
MatrixXf Lc, Bc;

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
			vPrimes.push_back(Vector3f(x, y, z));
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
vector<int> loadedVertTypes;

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

			loadedVertTypes.push_back(num);
		}

	}

	return true;
}

void initialRotations() {
	//For debugging
	initialRotation <<
		1, 0, 0,
		0, 1, 0,
		0, 0, 1;

	for (int i = 0; i < loadedVertices.size(); i++) {
		rotations.push_back(initialRotation);
	}
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

	/*printf("\nDeformation Matrix:\n");
	cout << deformationMatrix << endl;
	printf("-------------------------------------------\n");*/
	return true;
}

MatrixXf bMatrix;

void intializeBMatrix() {
	printf("Initializing B Matrix\n");
	Bc = MatrixXf(loadedVertices.size(), 3);
	Bc.fill(0);

	for (int r = 0; r < loadedVertTypes.size(); r++) {
		int vType = loadedVertTypes.at(r);
		if (vType == 0) {
			Bc.row(r) = loadedVertices.at(r);
		}
		else if (vType == 1) {
			//Free
		}
		else if (vType == 2) {
			Vector3f v3 = loadedVertices.at(r);

			Vector4f v4 = deformationMatrix * Vector4f(v3(0), v3(1), v3(2), 1);
			//cout << r << " " <<  v4 << endl;
			Bc.row(r) = Vector3f(v4(0), v4(1), v4(2));
		}
	}

	

	/*printf("\nBc Matrix:\n");
	cout << Bc << endl;
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

	for (int i = 0; i < numVertices; i++) {
		degMatrix(i, i) = neighbours[i].size();
	}

	/*printf("\nDegree Matrix (from adjacency):\n");
	cout << degMatrix << endl;
	printf("-------------------------------------------\n");*/

	laplacianMatrix = degMatrix.cast<float>() - adjMatrix.cast<float>();

	/*printf("\nLaplacian Matrix (L = D - A):\n");
	cout << laplacianMatrix << endl;
	printf("-------------------------------------------\n");*/

}

void modifyLaplacian() {
	printf("Modify Laplacian\n");

	//Constrain v wrt laplacian
	Lc = laplacianMatrix;


	for (int row = 0; row < loadedVertTypes.size(); row++) {
		int vertType = loadedVertTypes.at(row);
		
		if (vertType != 1) {
			for (int col = 0; col < loadedVertTypes.size(); col++) {
				Lc(row, col) = 0;
				if (row == col) {
					Lc(row, col) = 1;
				}
			}
		}
	}

	/*printf("\nLc Matrix:\n");
	cout << Lc << endl;
	printf("-------------------------------------------\n");*/
}

void updateBMatrix() {
	for (int i = 0; i < loadedVertices.size(); i++) {

		if (loadedVertTypes.at(i) == 1) {

			Bc.row(i) = Vector3f(0, 0, 0);
			for (int vNIndex : neighbours[i]) {

				Vector3f current = loadedVertices.at(i);

				Vector3f neighbour = loadedVertices.at(vNIndex);

				Vector3f newVec(.5 * (rotations.at(i) + rotations.at(vNIndex)) * (current - neighbour));
				//Check that this should be  += not =
				Bc.row(i) += newVec;
			}
		}
	}


	/*printf("\nBc Matrix Update:\n");
	cout << Bc << endl;
	printf("-------------------------------------------\n");*/
}

MatrixXf Dc;

void solvePositions() {
	printf("Solving Positions\n");
	Dc = Lc.inverse() * Bc;
	//Dc = Lc.householderQr().solve(Bc);
	
	/*printf("\nDc Matrix:\n");
	cout << Dc << endl;
	printf("-------------------------------------------\n");*/

	for (int i = 0; i < Dc.rows(); i++) {
		vPrimes.at(i) = Dc.row(i);
	}
}

void solveRotations() {
	printf("Solving Rotations\n");

	for (int i = 0; i < loadedVertices.size(); i++) {
		//printf("\n\n");
		set<int> nForI = neighbours[i];
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

		MatrixXf covariance = vHat * primeHat.transpose();
		JacobiSVD<MatrixXf> svd(covariance, ComputeThinU | ComputeThinV);
		MatrixXf uMatrix = svd.matrixU(), vMatrix = svd.matrixV();
		MatrixXf optimalRot = vMatrix * uMatrix.transpose();
		//cout << optimalRot << endl;
		rotations.at(i) = optimalRot;

		/*printf("\nRotation at %i:\n", i);
		cout << optimalRot << endl;
		printf("-------------------------------------------\n");*/
	}
}

//Exports the model after deformation
bool exportModel(int iteration) {
	printf("Exporting...\n");

	string path = "./ExportObj/deformedMesh" + to_string(iteration) + ".obj";

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
		//cout << "v " + to_string(v(0)) + " " + to_string(v(1)) + " " + to_string(v(2)) << endl;
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

void initialGuess() {

	int n = loadedVertTypes.size();
	int nc = 0;
	for (int row = 0; row < loadedVertTypes.size(); row++) {
		int vertType = loadedVertTypes.at(row);
		
		if (vertType != 1) {
			nc++;
		}

	}
	int nTotal = n + nc;	

	MatrixXf vInit(loadedVertices.size(), 3);
	for (int i = 0; i < loadedVertices.size(); i++) {
		vInit.row(i) = loadedVertices.at(i);
	}
	MatrixXf bInit = laplacianMatrix * vInit;

	MatrixXf LcInit(nTotal, nTotal);
	MatrixXf BcInit(nTotal, 3);

	//Constrain v wrt laplacian
	LcInit.block(0, 0, n, n) = laplacianMatrix;
	BcInit.block(0, 0, n, 3) = bInit;

	int cIx = 0;
	for (int vIx = 0; vIx < loadedVertTypes.size(); vIx++) {
		int vertType = loadedVertTypes.at(vIx);
		if (vertType == 0) {
			LcInit(n + cIx, vIx) = 1.0;
			LcInit(vIx, n + cIx) = 1.0;
			BcInit.row(n + cIx) = 1.0 * loadedVertices.at(vIx);
			cIx++;
		}
		else if (vertType == 2) {
			Vector3f v3 = loadedVertices.at(vIx);
			Vector4f v4d = deformationMatrix * Vector4f(v3(0), v3(1), v3(2), 1.0f);
			LcInit(n + cIx, vIx) = 1.0;
			LcInit(vIx, n + cIx) = 1.0;
			BcInit.row(n + cIx) = 1.0 * Vector3f(v4d(0), v4d(1), v4d(2));
			cIx++;
		}
	}

	MatrixXf vPositions = LcInit.inverse() * BcInit;

	for (int i = 0; i < vPrimes.size(); i++) {
		vPrimes.at(i) = vPositions.row(i);
	}

}

//Runs on startup
int main(int argc, char **argv) {

	//Load Files
	loadModel(meshFile);
	loadVertexTypes(vTypesFile);
	loadDeformation(deformationFile);

	generateLaplacian();

	modifyLaplacian();

	initialGuess();
	exportModel(-1);

	initialRotations();
	solveRotations();

	intializeBMatrix();

	for (int i = 0; i < iterations; i++) {
		printf("Iteration: %i\n", i);

		updateBMatrix();

		//solve for best positions
		solvePositions();

		//solve for best rotations
		solveRotations();

		//Export Model
		exportModel(i);
	}

	exit(0);
}