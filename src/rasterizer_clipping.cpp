#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::ivec2;
using glm::vec2;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

//Screen information
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;

//Camera information
vec3 cameraPos(0, 0, -3.001);
mat3 cameraRot(vec3(1,0,0), vec3(0,1,0), vec3(0,0,1));
float yaw = 0;
float pitch = 0;
float focalLength = 500;
float posDelta = 0.01;
float rotDelta = 0.01;

float clipBoundary = 0;
float maxDepth = 6.0f;

//Scene information
vector<Triangle> triangles;

//The depth buffer for the screen
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

//Floating point inaccuracy constant
float epsilon = 0.00001;

//Data structure holding information for each pixel
struct Pixel {
	int x;
	int y;
	float zinv;
	vec3 pos3d;
};

struct Vertex {
	vec3 o;	//position in original coordinate system
	vec3 c; //position in current coordinate system
	float w; 
};

vec3 lightPos(0,-0.5, -0.7);
float lightPosStep = 0.1f;
vec3 lightPower = 14.0f * vec3(1,1,1);
vec3 indirectLightPowerPerArea = 0.5f * vec3(1,1,1);

vec3 currentNormal;
vec3 currentReflectance;
/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void updateRotationMatrix();

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void VertexShader(const Vertex& v, Pixel& p);
void DrawPolygon(vector<Vertex> vertices);

void PixelShader(const Pixel& p);

int main( int argc, char* argv[] )
{
	//Load in the scene
	LoadTestModel(triangles);

	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );

	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	printf("Render time: %f ms.\n", dt);

	//get key presses and update camera position
	Uint8* keystate = SDL_GetKeyState(0);
	if(keystate[SDLK_UP])
	{
		//calculate the z-axis vector and move camera in the positive direction
		vec3 forward(cameraRot[2][0], cameraRot[2][1], cameraRot[2][2]);
		cameraPos += posDelta * forward;
	}
	if(keystate[SDLK_DOWN])
	{	
		//calculate the z-axis vector and move camera in the negative direction
		vec3 forward(cameraRot[2][0], cameraRot[2][1], cameraRot[2][2]);
		cameraPos -= posDelta * forward;
		
	}
	if(keystate[SDLK_LEFT])
	{
		//decrease the rotation angle and update rotation matrix 
		yaw -= rotDelta;
		updateRotationMatrix();
	}
	if(keystate[SDLK_RIGHT])
	{
		//increase the rotation angle and update the rotation matrix
		yaw += rotDelta;
		updateRotationMatrix();
	}
	if(keystate[SDLK_m])
	{
		pitch += rotDelta;
		updateRotationMatrix();
	}
	if(keystate[SDLK_n])
	{
		pitch -= rotDelta;
		updateRotationMatrix();
	}

	//Move light position depending on key press
	if(keystate[SDLK_w])
	{
		lightPos.z += lightPosStep;
	}
	if(keystate[SDLK_s])
	{
		lightPos.z -= lightPosStep;
	}
	if(keystate[SDLK_a])
	{
		lightPos.x -= lightPosStep;
	}
	if(keystate[SDLK_d])
	{
		lightPos.x += lightPosStep;
	}
	if(keystate[SDLK_q])
	{
		lightPos.y -= lightPosStep;
	}
	if(keystate[SDLK_e])
	{
		lightPos.y += lightPosStep;
	}
}

void Draw()
{
	SDL_FillRect(screen, 0, 0);

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	//reset the depth buffer

	for(int y = 0; y < SCREEN_HEIGHT; y++) {
		for(int x = 0; x < SCREEN_WIDTH; x++) {
			depthBuffer[y][x] = 0;
		}
	}

	//iterate through all of the triangles
	for(unsigned int i = 0; i < triangles.size(); ++i)
	{

		//Create a list of the triangles vertices
		vector<Vertex> vertices(3);
		//initalise the vertexes
		vertices[0].o = triangles[i].v0;
		vertices[1].o = triangles[i].v1;
		vertices[2].o = triangles[i].v2;
		
		currentNormal = triangles[i].normal;
		currentReflectance = triangles[i].color;

		//draw the triangle
		DrawPolygon(vertices);

	}


	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

//Update the cameras coordinate system
void updateRotationMatrix() {
	mat3 yRot(vec3(cos(yaw),0,-sin(yaw)), vec3(0,1,0), vec3(sin(yaw),0,cos(yaw)));
	mat3 xRot(vec3(1,0,0), vec3(0, cos(pitch), sin(pitch)), vec3(0, -sin(pitch), cos(pitch)));

	//Calcuate new columns for the camera's rotation matrix
	cameraRot = yRot * xRot;
}

int calculateOctant(Pixel a, Pixel b) {
	int dx = b.x - a.x;
	int dy = b.y - a.y;
	float gradient;
	if(dy == 0) {
		gradient = 0;
	}
	else if(dx == 0) {
		gradient = 10000;
	}
	else {
		gradient = (float) dy / (float) dx;
	}
	/*printf("(%d,%d), (%d,%d)\n", a.x, a.y, b.x, b.y);
	printf("gradient: %f\n", gradient);*/
	
	if(a.y < b.y && gradient > 1 && gradient < numeric_limits<int>::max()) {
		return 1;
	}
	else if(a.x < b.x && gradient >= 0 && gradient <= 1) {
		return 0;
	}
	else if(a.x < b.x && gradient < 0 && gradient >= -1) {
		return 7;
	}
	else if(a.y > b.y && gradient < -1 && gradient > numeric_limits<int>::min()) {
		return 6;
	}
	else if(a.y > b.y && gradient  > 1 && gradient  < numeric_limits<int>::max()) {
		return 5;
	}
	else if(a.x > b.x && gradient > 0 && gradient <= 1) {
		return 4;
	}
	else if(a.x > b.x && gradient <= 0 && gradient >= -1) {
		return 3;
	}
	else if(a.y < b.y && gradient < -1 && gradient > numeric_limits<int>::min()) {
		return 2;
	}
	else {
		return -1;
	}
}

void switchToOctantZeroFrom(int octant, Pixel& p) {
	int tmp = p.x;
	switch(octant) {
		case 0:
			break; 
		case 1:
			p.x = p.y;
			p.y = tmp;
			break;
		case 2:
			p.x = p.y;
			p.y = -tmp;
			break;
		case 3:
			p.x = -p.x;
			break;
		case 4:
			p.x = -p.x;
			p.y = -p.y;
			break;
		case 5:
			p.x = -p.y;
			p.y = -tmp;
			break;
		case 6:
			p.x = -p.y;
			p.y = tmp;
			break;
		case 7:
			p.y = -p.y;
			break;
	}
}

void switchFromOctantZeroTo(int octant, Pixel& p) {
	int tmp = p.x;
	switch(octant) {
		case 0:
			break;
		case 1:
			p.x = p.y;
			p.y = tmp;
			break;
		case 2:
			p.x = -p.y;
			p.y = tmp;
			break;
		case 3:
			p.x = -p.x;
			break;
		case 4:
			p.x = -p.x;
			p.y = -p.y;
			break;
		case 5:
			p.x = -p.y;
			p.y = -tmp;
			break;
		case 6:
			p.x = p.y;
			p.y = -tmp;
			break;
		case 7:
			p.y = -p.y;
			break;
	}
}

//Interpolate the points on a line between vectors a and b and put the points into result
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result) {

	int octant = calculateOctant(a, b);
	//printf("octant: %d\n", octant);
	switchToOctantZeroFrom(octant, a);
	switchToOctantZeroFrom(octant, b);

	int N = result.size();

	//printf("new coordinates: (%d,%d) (%d,%d)\n", a.x, a.y, b.x, b.y);
	float stepZinv = (b.zinv - a.zinv) / float(max(N-1,1));
	vec3 stepPos3d = (b.pos3d - a.pos3d) / float(max(N-1,1));

	float currentZinv = (float) a.zinv;
	vec3 currentPos3d = a.pos3d;

	//Calculate the steps needed in the x and y direction, and also the step needed for the zinv
	int dx = b.x - a.x;
	int dy = b.y - a.y;
	int D = 2 * dy - dx;
	int y = a.y;

	for(int i = 0; i < N; i++) {
		int x = i + a.x;
		result[i].x = x;
		result[i].y = y;
		result[i].zinv = currentZinv;
		result[i].pos3d = currentPos3d;
		if(D > 0) {
			y = y + 1;
			D = D - 2 * dx;
		}
		D = D + 2 * dy;
		currentZinv += stepZinv;
		currentPos3d += stepPos3d;
	}

	for(unsigned int i = 0; i < result.size(); i++) {
		switchFromOctantZeroTo(octant, result[i]);
	}
}

//Calculate the leftPixels and rightPixels lists from the vertex coordinats of the polygon
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
		
	//The maximum and minimum integer values
	int minY = numeric_limits<int>::max();
	int maxY = numeric_limits<int>::min();

	//Calculate the smallest and largest y values of the polygon
	for(unsigned int i = 0; i < vertexPixels.size(); i++) 
	{
		if(vertexPixels[i].y < minY)
			minY = vertexPixels[i].y;
		if(vertexPixels[i].y > maxY)
			maxY = vertexPixels[i].y;
	}

	//resize the leftPixels and rightPixels lists to the correct size
	unsigned int numRows = maxY - minY + 1;
	leftPixels.resize(numRows);
	rightPixels.resize(numRows);

	//set the leftPixels x values to the maximum integer value and
	//the rightPixels x values to the minimum integer value
	for(unsigned int i = 0; i < numRows; i++)
	{
		leftPixels[i].x = numeric_limits<int>::max();
		rightPixels[i].x = numeric_limits<int>::min();
	}


	//Iterate through all of the vertexes
	int V = vertexPixels.size();
	for(int i = 0; i < V; ++i)
	{
		//i is the index of the current vertex and j is the index of the next vertex
		int j = (i+1)%V;
		
		//a is the current vertex and b is the next vertex
		Pixel a = vertexPixels[i];
		Pixel b = vertexPixels[j];

		//calculate the number of pixels needed for the line between a and b		
		int pixels = max(abs(a.x - b.x), abs(a.y - b.y)) + 1;
		vector<Pixel> line(pixels);

		//Interpolate the points on the line between a and b
		Interpolate(a,b,line);
		/*for(unsigned int i = 0; i < line.size(); i++) {
			printf("(%d,%d) ", line[i].x, line[i].y);
		}
		printf("\n");*/
	
		//Iterate through all of the interpolated points on the line		
		for(unsigned int k = 0; k < line.size(); k++)
		{

			//if the x value for this y coordinate is less than the current x value for
			//this y coordinate then update the entry with the current point values
			if(leftPixels[line[k].y - minY].x > line[k].x)
			{
				leftPixels[line[k].y - minY].x = line[k].x;
				leftPixels[line[k].y - minY].y = line[k].y;
				leftPixels[line[k].y - minY].zinv = line[k].zinv;
				leftPixels[line[k].y - minY].pos3d = line[k].pos3d;
			}

			//if the x value for this y coordinate is greater than the current x value for
			//this y coordinate then update the entry with the current point values
			if(rightPixels[line[k].y - minY].x < line[k].x)
			{
				rightPixels[line[k].y - minY].x = line[k].x;
				rightPixels[line[k].y - minY].y = line[k].y;
				rightPixels[line[k].y - minY].zinv = line[k].zinv;
				rightPixels[line[k].y - minY].pos3d = line[k].pos3d;
			}
		}
	}
}

//Draw the polygon given the leftmost and rightmost pixel positions for each y value
void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels) {
	//Iterate through all the leftPixels and rightPixel elements
	unsigned int V = leftPixels.size();
	for(unsigned int i = 0; i < V; i++) {

		//Calculate the number of pixels needed for the line between corresponding leftPixels and rightPixels elements
		int pixels = rightPixels[i].x - leftPixels[i].x + 1;
		vector<Pixel> line(pixels);

		//Interpolate between these pixels (as we need to calculate the zinv values for each pixel)
		Interpolate(leftPixels[i], rightPixels[i], line);

		//Iterate over all the interpolated pixels
		for(int j = 0; j < pixels; j++) {
			//If the pixels zinv value is greater than the corresponding one in the buffer, then draw the pixel
			PixelShader(line[j]);
		}
	}
}

void DrawLine(Pixel a, Pixel b)
{
	int pixels = max(abs(a.x - b.x), abs(a.y - b.y)) + 1;
	vector<Pixel> line(pixels);
	Interpolate(a,b,line);
	for(unsigned int i = 0; i < line.size(); ++i) {
		PixelShader(line[i]);
	}
}

//Calculate the Euclidean distance between the two given vectors
float distanceBetweenPoints(vec3 a, vec3 b) {
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}

//Caluclate the dot product between the two given vectors
float dotProduct(vec3 a, vec3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

//Calculate the normalised direction vector from the given vector to the lightsource
vec3 unitVectorToLightSource(vec3 a) {
	vec3 v(lightPos.x - a.x, lightPos.y - a.y, lightPos.z - a.z);
	return normalize(v);
}

void TransformVertex(Vertex& v) {
	v.c = (v.o - cameraPos) * cameraRot;
}

//Project the 3d point v to 2D
void VertexShader(const Vertex& v, Pixel& p) {

	//project the x and y values
	p.x = round(focalLength * (v.c.x / v.c.z) + ((float) SCREEN_WIDTH / 2.0f));
	p.y = round(focalLength * (v.c.y / v.c.z) + ((float) SCREEN_HEIGHT / 2.0f));

	//calculate the inverse of the depth of the point
	p.zinv = 1.0f/(float)v.c.z;

	p.pos3d = v.o * p.zinv;
}

void PixelShader(const Pixel& p) {
	
	//Pi constant
	const float pi = 3.1415926535897;

	vec3 pos3d = p.pos3d / p.zinv;

	//distance from intersection point to light source
	float radius = distanceBetweenPoints(pos3d, lightPos);

	//The power per area at this point
	vec3 B = lightPower / (4 * pi * (float) pow(radius,3));

	//unit vector describing normal of surface
	vec3 n = currentNormal;

	//unit vector describing direction from surface point to light source
	vec3 r = unitVectorToLightSource(pos3d);

	//fraction of the power per area depending on surface's angle from light source
	vec3 D = B * max(dotProduct(r,n),0.0f);

	vec3 illumination = currentReflectance * (D + indirectLightPowerPerArea);

	//If pixels depth is less than the current pixels depth in the image
	//then update the image
	if(p.x >= 0 && p.x < SCREEN_WIDTH && p.y >= 0 && p.y < SCREEN_HEIGHT) {
		if(p.zinv > depthBuffer[p.y][p.x] + epsilon && p.zinv > 0) {
			depthBuffer[p.y][p.x] = p.zinv;
			PutPixelSDL(screen, p.x, p.y, illumination);
		}
	}
}

bool Intersection(const Vertex a, const Vertex b) {
	float divide = b.c.x - a.c.x;
	if(divide < epsilon && divide > -epsilon) {
		return false;
	}
	else {
		return true;
	}
}

void ClipSpace(vector<Vertex>& vertices) {
	for(unsigned int i = 0; i < vertices.size(); i++) {
		vertices[i].w = vertices[i].c.z / focalLength;
	}
}

Vertex ClipRight(Vertex start, Vertex end) {
	float factor = (((float) SCREEN_WIDTH) / 2.0f) - clipBoundary;

	float a = (start.c.x - start.w * factor) / (-start.w * factor + end.w * factor + start.c.x - end.c.x);
	Vertex P;
	P.c = (1.0f - a) * start.c + a * end.c;
	P.o = (1.0f - a) * start.o + a * end.o;

	float w = (1.0f - a) * start.w + a * end.w;
	float ratio = w / (P.c.z / focalLength);
	P.c = P.c / ratio;
	P.w = (P.c.z / focalLength);

	return P;
}

void ClipRightEdge(vector<Vertex>& inputList, vector<Vertex>& outputList) {

    Vertex start = inputList[inputList.size() - 1];
	for(unsigned int i = 0; i < inputList.size(); i++) {
		Vertex end = inputList[i];

		float startxmax = start.w * ((float) SCREEN_WIDTH / 2.0f - clipBoundary);
		float endxmax = end.w * ((float) SCREEN_WIDTH / 2.0f - clipBoundary);

		if(end.c.x < endxmax) {
			if(start.c.x > startxmax) {
				Vertex P = ClipRight(start, end);
				outputList.push_back(P);
			}
			outputList.push_back(end);
		}
		else if(start.c.x < startxmax) {
			Vertex P = ClipRight(start, end);
			outputList.push_back(P);
		}
		start = end;
	}
}

Vertex ClipLeft(Vertex start, Vertex end) {
	float factor = (((float) SCREEN_WIDTH) / 2.0f) - clipBoundary;

	float a = (start.w * factor + start.c.x) / ((start.w * factor + start.c.x) - (end.w * factor + end.c.x));
	Vertex P;
	P.c = (1.0f - a) * start.c + a * end.c;
	P.o = (1.0f - a) * start.o + a * end.o;

	float w = (1.0f - a) * start.w + a * end.w;
	float ratio = w / (P.c.z / focalLength);
	P.c = P.c / ratio;
	P.w = (P.c.z / focalLength);

	return P;
}

void ClipLeftEdge(vector<Vertex>& inputList, vector<Vertex>& outputList) {

    Vertex start = inputList[inputList.size() - 1];
	for(unsigned int i = 0; i < inputList.size(); i++) {
		Vertex end = inputList[i];

		float startxmin = -start.w * ((float) SCREEN_WIDTH / 2.0f - clipBoundary);
		float endxmin = -end.w * ((float) SCREEN_WIDTH / 2.0f - clipBoundary);

		if(end.c.x > endxmin) {
			if(start.c.x < startxmin) {
				Vertex P = ClipLeft(start, end);
				outputList.push_back(P);
			}
			outputList.push_back(end);
		}
		else if(start.c.x > startxmin) {
			Vertex P = ClipLeft(start, end);
			outputList.push_back(P);
		}
		start = end;
	}
}

Vertex ClipTop(Vertex start, Vertex end) {
	float factor = (((float) SCREEN_HEIGHT) / 2.0f) - clipBoundary;

	float a = (start.w * factor + start.c.y) / ((start.w * factor + start.c.y) - (end.w * factor + end.c.y));
	Vertex P;
	P.c = (1.0f - a) * start.c + a * end.c;
	P.o = (1.0f - a) * start.o + a * end.o;

	float w = (1.0f - a) * start.w + a * end.w;
	float ratio = w / (P.c.z / focalLength);
	P.c = P.c / ratio;
	P.w = (P.c.z / focalLength);

	return P;
}

void ClipTopEdge(vector<Vertex>& inputList, vector<Vertex>& outputList) {

    Vertex start = inputList[inputList.size() - 1];
	for(unsigned int i = 0; i < inputList.size(); i++) {
		Vertex end = inputList[i];

		float startymin = -start.w * ((float) SCREEN_HEIGHT / 2.0f - clipBoundary);
		float endymin = -end.w * ((float) SCREEN_HEIGHT / 2.0f - clipBoundary);

		if(end.c.y > endymin) {
			if(start.c.y < startymin) {
				Vertex P = ClipTop(start, end);
				outputList.push_back(P);
			}
			outputList.push_back(end);
		}
		else if(start.c.y > startymin) {
			Vertex P = ClipTop(start, end);
			outputList.push_back(P);
		}
		start = end;
	}
}

Vertex ClipBottom(Vertex start, Vertex end) {
	float factor = (((float) SCREEN_HEIGHT) / 2.0f) - clipBoundary;

	float a = (start.c.y - start.w * factor) / (-start.w * factor + end.w * factor + start.c.y - end.c.y);
	Vertex P;
	P.c = (1.0f - a) * start.c + a * end.c;
	P.o = (1.0f - a) * start.o + a * end.o;

	float w = (1.0f - a) * start.w + a * end.w;
	float ratio = w / (P.c.z / focalLength);
	P.c = P.c / ratio;
	P.w = (P.c.z / focalLength);

	return P;
}

void ClipBottomEdge(vector<Vertex>& inputList, vector<Vertex>& outputList) {

    //Vertex start = inputList[inputList.size() - 1];
	for(unsigned int i = 0; i < inputList.size(); i++) {
		Vertex start = inputList[i];
		Vertex end = inputList[(i+1)%inputList.size()];

		float startymax = start.w * ((float) SCREEN_HEIGHT / 2.0f - clipBoundary);
		float endymax = end.w * ((float) SCREEN_HEIGHT / 2.0f - clipBoundary);

		if(end.c.y < endymax) {
			if(start.c.y > startymax) {
				Vertex P = ClipBottom(start, end);
				outputList.push_back(P);
			}
			outputList.push_back(end);
		}
		else if(start.c.y < startymax) {
			Vertex P = ClipBottom(start, end);
			outputList.push_back(P);
		}
		//start = end;
	}
}


bool Clip(vector<Vertex>& vertices) {

	ClipSpace(vertices);

	vector<Vertex> outputList = vertices;
	vector<Vertex> inputList;

	inputList = outputList;
	outputList.clear();

	ClipBottomEdge(inputList, outputList);

	inputList = outputList;
	outputList.clear();

	ClipRightEdge(inputList, outputList);

	inputList = outputList;
	outputList.clear();

	ClipTopEdge(inputList, outputList);

	inputList = outputList;
	outputList.clear();

	ClipLeftEdge(inputList, outputList);

	vertices = outputList;

	return (vertices.size() > 0);
}

bool IsPolygonInfrontOfCamera(const vector<Vertex>& vertices) {
	bool visible = false;
	for(unsigned int i = 0; i < vertices.size(); i++) {
		if(vertices[i].c.z > 0) {
			visible = true;
		}
	}
	return visible;
}

bool IsPolygonWithinMaxDepth(const vector<Pixel>& vertexPixels) {
	float invMaxDepth = 1.0f / maxDepth;

	bool visible = false;
	for(unsigned int i = 0; i < vertexPixels.size(); i++) {
		if(vertexPixels[i].zinv > invMaxDepth) {
			visible = true;
		}
	}
	return visible;
}

//Draw a polygon given its vertices, taking into account depth
void DrawPolygon(vector<Vertex> vertices)
{

	//Transform world
	for(unsigned int i = 0; i < vertices.size(); i++) {
		TransformVertex(vertices[i]);
	}

	//Backface culling
	if(IsPolygonInfrontOfCamera(vertices)) {

		//Clip polygon - function returns false if new polygon has zero vertices
		if(Clip(vertices)) {

			//Calculate the projection of the polygons vertexes
			vector<Pixel> vertexPixels(vertices.size());
			for(unsigned int i = 0; i < vertices.size(); i++) {
				VertexShader(vertices[i], vertexPixels[i]);
			}

			if(IsPolygonWithinMaxDepth(vertexPixels)) {
		
				//lists to store the left most and right post pixel x values for each y value
				vector<Pixel> leftPixels;
				vector<Pixel> rightPixels;

				//Compute the leftPixels and rightPixels lists from the vertexes
				ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
				//printf("compute polygon rows done\n");

				//Draw the polygon
				DrawPolygonRows(leftPixels, rightPixels);
				//printf("draw polygon rows done\n");
			}
		}
	}
}