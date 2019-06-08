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

//==========================================================================================
// GLOBAL VARIABLES

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
const float focalLength = 500;
const float posDelta = 0.02;
const float rotDelta = 0.02;

//Lighting and colour information
vec3 lightPos(0,-0.5, -0.7);
const float lightPosStep = 0.1f;
const vec3 lightPower = 14.0f * vec3(1,1,1);
const vec3 indirectLightPowerPerArea = 0.5f * vec3(1,1,1);
vec3 currentNormal;
vec3 currentReflectance;

//The depth buffer for the screen
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

//Floating point inaccuracy constant
const float epsilon = 0.00001;

//clipping information
const float clipBoundary = 5;
const float maxDepth = 100;

//Texture information
SDL_Surface* checkered512x512;
enum Texture {None, Checkered};
Texture currentTexture = None;

//rasterizer features
const bool clipping = true;
const bool textures = true;

//===============================================================================================
// DATA STRUCTURES

//Data structure holding information for each pixel
struct Pixel {
	int x;
	int y;
	float zinv;
	vec3 pos3d;
	ivec2 textureCoordinates;
};

struct Vertex {
	vec3 o;	//position in original coordinate system
	vec3 c; //position in current coordinate system
	float w;
	ivec2 textureCoordinates;
};

//==============================================================================================
// FUNCTIONS
void LoadTextures();
void Update();
void Draw(const vector<Triangle> triangles);
void updateRotationMatrix();
void Interpolate(const Pixel a, const Pixel b, vector<Pixel>& result);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void VertexShader(const Vertex& v, Pixel& p);
void DrawPolygon(const vector<Vertex> vertices);
void PixelShader(const Pixel p);


int main( int argc, char* argv[] )
{
	if(textures) {
	  //Load in the textures
	  LoadTextures();
	}

	//Scene information
	vector<Triangle> triangles;
	//Load in the scene
	LoadTestModel(triangles);

	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw(triangles);
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );

	return 0;
}

//Load the textures to use on surfaces
void LoadTextures() {
	checkered512x512 = SDL_LoadBMP("bin/images/checkered512x512.bmp");
}

//Update the camera and light positions
void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

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

//Rendering the screen for the current frame
void Draw(const vector<Triangle> triangles)
{
	//clear the screen
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
		//set the vertices original positions equal to the triangles vertices positions
		vertices[0].o = triangles[i].v0.pos3d;
		vertices[1].o = triangles[i].v1.pos3d;
		vertices[2].o = triangles[i].v2.pos3d;
		
		//set the current normal vector and reflectance colour equal to the triangles normal and colour respectively
		currentNormal = triangles[i].normal;
		currentReflectance = triangles[i].color;

		//If we are using textures then set the current texture and give the vertices their texture coordinates
		if(textures) {
			if(triangles[i].texture == 0) {
				currentTexture = None;
			} 
			else {
				if(triangles[i].texture == 1) {
					currentTexture = Checkered;
				}
				vertices[0].textureCoordinates = triangles[i].v0.textureCoordinates;
				vertices[1].textureCoordinates = triangles[i].v1.textureCoordinates;
				vertices[2].textureCoordinates = triangles[i].v2.textureCoordinates;
			}
		}
		
		//draw the triangle
		DrawPolygon(vertices);
	}


	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

    //update the screen
	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

//Update the cameras coordinate system
void updateRotationMatrix() {
	//create the rotation matrix for the rotation around the y axis
	mat3 yRot(vec3(cos(yaw),0,-sin(yaw)), vec3(0,1,0), vec3(sin(yaw),0,cos(yaw)));
	//create the rotation matrix for the rotation around the x axis
	mat3 xRot(vec3(1,0,0), vec3(0, cos(pitch), sin(pitch)), vec3(0, -sin(pitch), cos(pitch)));

	//Calcuate the new camera rotation matrix
	cameraRot = yRot * xRot;
}

//Interpolate the points on a line between vectors a and b and put the points into the result vector
void Interpolate(const Pixel a, const Pixel b, vector<Pixel>& result) {

	//Calculate the steps needed to interpolate between pixels
	int N = result.size();
	float stepX = (float) (b.x - a.x) / float(max(N-1,1));
	float stepY = (float) (b.y - a.y) / float(max(N-1,1));
	float stepZinv = (b.zinv - a.zinv) / float(max(N-1,1));
	vec3 stepPos3d = (b.pos3d - a.pos3d) / float(max(N-1,1));

	//Set the accumulator values to the same as vector a
	float currentX = (float) a.x;
	float currentY = (float) a.y;
	float currentZinv = (float) a.zinv;
	vec3 currentPos3d = a.pos3d;

	//Calculate each point on the line
	for(int i = 0; i < N; i++)
	{
		//Store the interpolated point in the result list
		result[i].x = round(currentX);
		result[i].y = round(currentY);
		result[i].zinv = currentZinv;
		result[i].pos3d = currentPos3d;

		//Update the accumulator values
		currentX += stepX;
		currentY += stepY;
		currentZinv += stepZinv;
		currentPos3d += stepPos3d;
	}

	//if textures are enabled
	if(textures) {
		//Interpolate the texture coordinates
		if(currentTexture != None) {
			for(int i = 0; i < N; i++) {
				float q = (float) i / (float) N;
				result[i].textureCoordinates.x = round( ( ( (float) a.textureCoordinates.x * a.zinv ) * (1 - q) + ( (float) b.textureCoordinates.x * b.zinv ) * q ) / result[i].zinv );
				result[i].textureCoordinates.y = round( ( ( (float) a.textureCoordinates.y * a.zinv ) * (1 - q) + ( (float) b.textureCoordinates.y * b.zinv ) * q ) / result[i].zinv );			
			}
		}
	}
}

//Calculate the leftPixels and rightPixels lists from the vertex coordinates of the polygon
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
				if(textures) {
					leftPixels[line[k].y - minY].textureCoordinates = line[k].textureCoordinates;
				}
			}

			//if the x value for this y coordinate is greater than the current x value for
			//this y coordinate then update the entry with the current point values
			if(rightPixels[line[k].y - minY].x < line[k].x)
			{
				rightPixels[line[k].y - minY].x = line[k].x;
				rightPixels[line[k].y - minY].y = line[k].y;
				rightPixels[line[k].y - minY].zinv = line[k].zinv;
				rightPixels[line[k].y - minY].pos3d = line[k].pos3d;
				if(textures) {
					rightPixels[line[k].y - minY].textureCoordinates = line[k].textureCoordinates;
				}
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

	p.textureCoordinates = v.textureCoordinates;
}

//calculate the colour needed for a pixel, draw it to the buffer and add its depth to the depth buffer
void PixelShader(const Pixel p) {
	
	//Pi constant
	const float pi = 3.1415926535897;

	vec3 pos3d = p.pos3d / p.zinv;

	//distance from intersection point to light source
	float radius = distanceBetweenPoints(pos3d, lightPos);

	//The power per area at this point
	vec3 B = lightPower / ((float) 4.0 * pi * (float) pow(radius,3));

	//unit vector describing normal of surface
	vec3 n = currentNormal;

	//unit vector describing direction from surface point to light source
	vec3 r = unitVectorToLightSource(pos3d);

	//fraction of the power per area depending on surface's angle from light source
	vec3 D = B * max(dotProduct(r,n),0.0f);

    //set the colour to the colour of current triangle
	vec3 color = currentReflectance;

    //if textures are enabled then set the colour equal to the colour at the texture coordinate
	if(textures) {
		if(currentTexture == Checkered) {
			if(p.textureCoordinates.x >= 0 && p.textureCoordinates.x <= 511 && p.textureCoordinates.y >= 0 && p.textureCoordinates.y <= 511) {
				color = GetPixelSDL(checkered512x512, p.textureCoordinates.x, p.textureCoordinates.y);
			}
		}
	}

    //change the colour due to light sources
	vec3 illumination = color * (D + indirectLightPowerPerArea);

	//If pixels depth is less than the current pixels depth in the image
	//then update the image
	if(p.x >= 0 && p.x < SCREEN_WIDTH && p.y >= 0 && p.y < SCREEN_HEIGHT) {
		if(p.zinv > depthBuffer[p.y][p.x] + epsilon && p.zinv > 0) {
			depthBuffer[p.y][p.x] = p.zinv;
			PutPixelSDL(screen, p.x, p.y, illumination);
		}
	}
}

//move all vertices into clip space
vector<Vertex> ClipSpace(const vector<Vertex> vertices) {
	vector<Vertex> newVertices = vertices;
	for(unsigned int i = 0; i < newVertices.size(); i++) {
		newVertices[i].w = newVertices[i].c.z / focalLength;
	}
	return newVertices;
}

Vertex ClipRight(Vertex start, Vertex end) {
	float factor = (((float) SCREEN_WIDTH) / 2.0f) - clipBoundary;

    //The percentage of how far from start to end the intersection ocurrs
	float a = (start.c.x - start.w * factor) / (-start.w * factor + end.w * factor + start.c.x - end.c.x);
	Vertex P;
	//The coordinates are a combination of the start and end coordinates
	P.c = (1.0f - a) * start.c + a * end.c;
	P.o = (1.0f - a) * start.o + a * end.o;
	if(textures) {
		P.textureCoordinates.x = (1.0f - a) * start.textureCoordinates.x + a * end.textureCoordinates.x;
		P.textureCoordinates.y = (1.0f - a) * start.textureCoordinates.y + a * end.textureCoordinates.y;
	}

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

		if(end.c.x < endxmax + epsilon) {
			if(start.c.x > startxmax - epsilon) {
				Vertex P = ClipRight(start, end);
				outputList.push_back(P);
			}
			outputList.push_back(end);
		}
		else if(start.c.x < startxmax + epsilon) {
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
	if(textures) {
		P.textureCoordinates.x = (1.0f - a) * start.textureCoordinates.x + a * end.textureCoordinates.x;
		P.textureCoordinates.y = (1.0f - a) * start.textureCoordinates.y + a * end.textureCoordinates.y;
	}

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

		if(end.c.x > endxmin - epsilon) {
			if(start.c.x < startxmin + epsilon) {
				Vertex P = ClipLeft(start, end);
				outputList.push_back(P);
			}
			outputList.push_back(end);
		}
		else if(start.c.x > startxmin - epsilon) {
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
	if(textures) {
		P.textureCoordinates.x = (1.0f - a) * start.textureCoordinates.x + a * end.textureCoordinates.x;
		P.textureCoordinates.y = (1.0f - a) * start.textureCoordinates.y + a * end.textureCoordinates.y;
	}

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

		if(end.c.y > endymin - epsilon) {
			if(start.c.y < startymin + epsilon) {
				Vertex P = ClipTop(start, end);
				outputList.push_back(P);
			}
			outputList.push_back(end);
		}
		else if(start.c.y > startymin - epsilon) {
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
	if(textures) {
		P.textureCoordinates.x = (1.0f - a) * start.textureCoordinates.x + a * end.textureCoordinates.x;
		P.textureCoordinates.y = (1.0f - a) * start.textureCoordinates.y + a * end.textureCoordinates.y;
	}

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

		if(end.c.y < endymax + epsilon) {
			if(start.c.y > startymax - epsilon) {
				Vertex P = ClipBottom(start, end);
				outputList.push_back(P);
			}
			outputList.push_back(end);
		}
		else if(start.c.y < startymax + epsilon) {
			Vertex P = ClipBottom(start, end);
			outputList.push_back(P);
		}
	}
}


vector<Vertex> Clip(const vector<Vertex> vertices) {

	vector<Vertex> clipSpaceVertices = ClipSpace(vertices);

	vector<Vertex> outputList = clipSpaceVertices;
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

	return outputList;
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
void DrawPolygon(const vector<Vertex> vertices)
{

	vector<Vertex> verticesCopy = vertices;

	//Transform world
	for(unsigned int i = 0; i < vertices.size(); i++) {
		TransformVertex(verticesCopy[i]);
	}

	//Backface culling
	if(IsPolygonInfrontOfCamera(verticesCopy)) {

		//If clipping is enabled then clip the polygon
		if(clipping) {
			verticesCopy = Clip(verticesCopy);
		}

		//Calculate the projection of the polygons vertexes
		vector<Pixel> vertexPixels(verticesCopy.size());
		for(unsigned int i = 0; i < verticesCopy.size(); i++) {
			VertexShader(verticesCopy[i], vertexPixels[i]);
		}

		if(IsPolygonWithinMaxDepth(vertexPixels)) {
		
			//lists to store the left most and right post pixel x values for each y value
			vector<Pixel> leftPixels;
			vector<Pixel> rightPixels;

			//Compute the leftPixels and rightPixels lists from the vertexes
			ComputePolygonRows(vertexPixels, leftPixels, rightPixels);

			//Draw the polygon
			DrawPolygonRows(leftPixels, rightPixels);
		}
	}
}